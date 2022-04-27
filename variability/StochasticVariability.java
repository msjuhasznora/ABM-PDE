package stochasticVariability;

// a refactored version of Sadegh's code

// We introduce a new cell type: a cell of "type 3" represents an infected cell that does not produce virus.
// The category is introduced in order to track the number of infected cells that got infected because of
// one specific cell.

// init description todo

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;
import HAL.Rand;

import java.io.File;
import java.util.Random;

import static HAL.Util.*;

public class StochasticVariability {

	public static int numberOfExperiments = 2000;
	public static int numberOfTicks = 7000; // per experiment

	// this data set will store the number of new infections at each time step for each experiment
	public static double[][] branchingProcessData = new double[numberOfExperiments][numberOfTicks];

	public static int x = 200;
	public static int y = 200;

	public static void main(String[] args) {

		java.util.Date now = new java.util.Date();
		java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		String date_time = dateFormat.format(now);

		String projPath = PWD() + "/stochasticVariability";
		final String output_dir = projPath + "/output/" + date_time + "/";

		new File(output_dir).mkdirs();
		FileIO outfile = new FileIO(output_dir.concat("/").concat("Out").concat(".csv"), "w");
		FileIO dataLastTick = new FileIO(output_dir.concat("/").concat("dataLastTick").concat(".csv"), "w");

		for (int experimentIterator = 0; experimentIterator < numberOfExperiments; experimentIterator++){

			System.out.println(experimentIterator);

			StochasticExperiment experiment = new StochasticExperiment(x, y);
			experiment.Init();

			double numberOfType3Cells = experiment.RunExperiment(numberOfTicks, experimentIterator);
			System.out.println(numberOfType3Cells);
		}

		for (int experimentIterator = 0; experimentIterator < branchingProcessData.length; experimentIterator++) {

			for (int tick = 0; tick < branchingProcessData[0].length ; tick++){
				outfile.Write(branchingProcessData[experimentIterator][tick] +",");
			}

			outfile.Write("\n");
			dataLastTick.Write(branchingProcessData[experimentIterator][branchingProcessData[0].length-1] + "\n");

		}

		outfile.Close();
		dataLastTick.Close();
	}
}


class StochasticExperiment extends AgentGrid2D<Cells>{

	public PDEGrid2D virusCon;
	public Rand rn;

	public double[] cellularVirusCon = new double[length];
	public double propHealthy = 0.9995;

	// SARS-CoV-2 parameters
	public double virusRemovalRate = 1.67 * Math.pow(10, -3);
	public double virusMax = 3.72 * Math.pow(10, -3);  // f_{i,j}
	public double virusDiffCoeff = 0.2; // D_V [sigma^2 / min]
	public double deathProb = 7.2 * Math.pow(10,-4);
	public double infectionRate = 1.01 * Math.pow(10,-7);

	public StochasticExperiment(int x, int y){
		super(x, y, Cells.class);
		this.rn = rn;
		virusCon = new PDEGrid2D(xDim, yDim);
		virusCon.Update();
	}

	void Init(){

		double randv = randomGenerator();

		double R = Math.round(randv * 10000.0);
		R = R/10000.0;

		for (int i = 0; i < length; i++){

			if( i == Math.round(R * xDim * yDim)) {
				Cells c = NewAgentSQ(i);
				c.CellInit(false,true,false,false);
			} else {
				Cells c = NewAgentSQ(i);
				c.CellInit(true,false,false,false);
			}
		}
	}

	double RunExperiment(int numberOfTicks, int experimentIterator){

		double[] cellCount = CountCells();

		for (int tick = 0; tick < numberOfTicks; tick ++){

			TimeStep();

			double totalVirusCon = TotalVirusCon();
			cellCount = CountCells();

			StochasticVariability.branchingProcessData[experimentIterator][tick] = cellCount[3];

			if (totalVirusCon < 0.001) {

				for (int t = tick; t < StochasticVariability.branchingProcessData[0].length-1 ; t++){
					StochasticVariability.branchingProcessData[experimentIterator][t+1] = StochasticVariability.branchingProcessData[experimentIterator][t];
				}

				tick = numberOfTicks;
			}
		}

		return cellCount[3];
	}

	double[] CountCells(){

		double healthyCells = 0, infectedCells = 0, deadCells = 0, type3Cells = 0;
		double[] cellCount = new double[4];

		for (Cells cell: this){
			if (cell.CellType == 0){
				healthyCells += 1;
			} else if (cell.CellType == 1 ){
				infectedCells += 1;
			} else if (cell.CellType == 2){
				deadCells += 1;
			} else if (cell.CellType == 3){
				type3Cells += 1;
			}

		}

		cellCount[0] = healthyCells;
		cellCount[1] = infectedCells;
		cellCount[2] = deadCells;
		cellCount[3] = type3Cells;

		return cellCount;

	}

	public Random rnd = new Random();
	double randomGenerator() {
		return rnd.nextDouble();
	}

	public void TimeStep(){

		TimeStepVirus();
		TimeStepCells();

	}

	void TimeStepCells(){

		for (Cells cell : this){
			cell.CellStepInfection();
		}

		for (Cells cell : this) {
			cell.CellStepDeath();
		}

	}

	void TimeStepVirus(){

		// decay of the virus
		for (Cells cell : this){
			virusCon.Add(cell.Isq(), -virusRemovalRate * virusCon.Get(cell.Isq()));
		}
		virusCon.Update();

		// virus production
		for (Cells cell : this){
			if (cell.CellType == 1){ // infected cell
				double addedVirusCon = VirusSource();
				double currentVirusCon = virusCon.Get(cell.Isq());
				double newVirusCon = addedVirusCon + currentVirusCon;
				virusCon.Set(cell.Isq(), newVirusCon);
			}
		}

		virusCon.DiffusionADI(virusDiffCoeff);
		virusCon.Update();

	}

	double TotalVirusCon() {

		double totalVirusCon = 0;
		for (int i = 0; i < length; i++){
			cellularVirusCon[i] = virusCon.Get(i);
		}

		for (double virusConInCell : cellularVirusCon ){
			totalVirusCon = totalVirusCon + virusConInCell;
		}

		return totalVirusCon;

	}

	double VirusSource(){

		return virusMax;

	}

	public void DrawModel(GridWindow vis){

		for (int i = 0; i < length; i++) {

			Cells drawMe = GetAgent(i);

			if (drawMe == null){
				vis.SetPix(i, RGB256(255, 255, 255));
			} else{

				if (drawMe.CellType == 0){ // healthy cells
					vis.SetPix(i, RGB256(0, 255, 51));
				}
				else if (drawMe.CellType == 1){ // infected cells
					vis.SetPix(i, RGB256(255, 0, 0));
				}
				else if (drawMe.CellType == 2){
					vis.SetPix(i, RGB256(0, 0, 0));
				}
				else if (drawMe.CellType == 3){
					vis.SetPix(i, RGB256(180, 3, 255));
				}

				vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapRGB(virusCon.Get(i)));

			}
		}
	}

}

class Cells extends AgentSQ2Dunstackable<StochasticExperiment>{

	int CellType;

	// cell types: 0 - healthy, 1 - infected, 2 - dead, 3 - an infected cell that does not produce virus

	public void CellInit(boolean isHealthy, boolean isInfected, boolean isDead, boolean isType3){

		if(isHealthy == true){
			this.CellType = 0;
		} else if(isInfected == true){
			this.CellType = 1;
		} else if(isDead == true) {
			this.CellType = 2;
		} else if(isType3 == true) {
			this.CellType = 3;
		}
	}

	public void CellStepInfection(){

		double virusConAtCell = G.virusCon.Get(Isq());

		double infectionProb = G.infectionRate * G.xDim * G.yDim * virusConAtCell;
		if (G.randomGenerator() < infectionProb) {

			if (this.CellType == 0){ // healthy cell
				this.CellType = 3; // a "type 3" cell is an infected cell that does not produce virus
			}
		}
	}

	public void CellStepDeath(){

		if (this.CellType == 1) {

			if(G.randomGenerator() < G.deathProb){
				this.CellType = 2;
			}
		}
	}

}