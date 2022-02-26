package NirmatrelvirExperiments;

// Experiment nr 1: nirmatrelvir alone (MOA: reduced virus production)
// Experiment nr 2: ritonavir-boosted nirmatrelvir (MOA: reduced virus production combined with
// a slower metabolism / decay of the antiviral drug)

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Gui.GifMaker;
import HAL.Tools.FileIO;
import HAL.Rand;
import static HAL.Util.*;

import java.io.File;

public class NirmatrelvirExperiments extends AgentGrid2D<Cells>{
	public PDEGrid2D virusConcentration;
	public PDEGrid2D drugConcentration;
	public Rand rn;
	public double[] cellularVirusConcentration = new double[length];
	public double[] cellularDrugConcentration = new double[length];

	public double ratioHealthy = 0.998, ratioInfected = 0.002,
			ratioCapillaries = 0;
	public double virusRemovalRate = 0.15, virusMax = 0.17,
			drugMax = 0.3;
	public double infectionRate = 0.15 * Math.pow(10,-5);
	public double deathProb = 0.011;
	public double diffCoeff = 0.6;
	public boolean isRitonavirBoosted = true;
	public double drugDecay = 0.05;

	public double MAX_PDE_STEP = 1;
	public double threshold = 0.000001;

	public NirmatrelvirExperiments(int xDim, int yDim, Rand rn, boolean isRitonavirBoosted){
		super(xDim, yDim, Cells.class);
		this.rn = rn;
		this.isRitonavirBoosted = isRitonavirBoosted;
		virusConcentration = new PDEGrid2D(xDim, yDim);
		drugConcentration = new PDEGrid2D(xDim, yDim);
		virusConcentration.Update();
		drugConcentration.Update();
	}


	public static void main(String[] args) {
		int y = 200, x = 200, visScale = 2;
		int numberOfTicks = 3200;
		boolean isRitonavirBoosted = false;

		java.util.Date now = new java.util.Date();
		java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		String date_time = dateFormat.format(now);
		String projPath = PWD() + "/NirmatrelvirExperiments";
		String outputDir = projPath + "/output/" + date_time + "/";
		new File(outputDir).mkdirs();

		GridWindow win = new GridWindow("Cellular state space, virus concentration, drug concentration", x*3, y, visScale,true);
		FileIO outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
		FileIO paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");

		NirmatrelvirExperiments model = new NirmatrelvirExperiments(x, y, new Rand(1), isRitonavirBoosted);
		GifMaker gm = new GifMaker(outputDir.concat("/").concat("test").concat(".gif"),100,true);
		// win.ToGIF("test.jpg");

		paramFile.Write("Parameters: \n ratioHealthy, virusRemovalRate, diffCoeff \n");
		paramFile.Write(model.ratioHealthy + "," + model.virusRemovalRate + "," + model.diffCoeff + "\n");
		outFile.Write("tick, healthyCells, infectedCells, deadCells, "
				+ "integratedVirusLoad, integratedDrugLoad \n");
		model.Init(); // Put X for defining boundaries
		model.DrawModel(win);
		gm.AddFrame(win);

		double[] cellCounts = model.CountCells();
		System.out.println(cellCounts[0]+", " + cellCounts[1] + ", " + cellCounts[2]);

		for (int tick = 0; tick < numberOfTicks; tick ++){
			System.out.println(tick);
			model.ModelStep(tick);
			model.DrawModel(win);
			// gm.AddFrame(win);
			double totalVirusConcentration = model.TotalVirusConcentration();
			double totalDrugConcentration = model.TotalDrugConcentration();
			cellCounts = model.CountCells();
			outFile.Write(tick +"," + cellCounts[0] + "," + cellCounts[1]+
					"," + cellCounts[2]+"," + cellCounts[3] + "," + totalVirusConcentration + "," + totalDrugConcentration + "\n");

		}

		gm.Close();
		outFile.Close();
		paramFile.Close();
		win.Close();
	}

	double[] CountCells(){
		double healthyCells = 0, infectedCells = 0, deadCells = 0,
				capillaryCells = 0;
		double[] cellCount = new double[4];
		for (Cells cell: this){
			if (cell.CellType == 0){
				healthyCells += 1;
			}
			else if (cell.CellType == 1 ){
				infectedCells += 1;
			}
			else if (cell.CellType == 2){
				deadCells += 1;
			}
			else if (cell.CellType == 3){
				capillaryCells += 1;
			}
		}

		cellCount[0] = healthyCells;
		cellCount[1] = infectedCells;
		cellCount[2] = deadCells;
		cellCount[3] = capillaryCells;

		return cellCount;
	}

	public void Init(){

		if (isRitonavirBoosted == true){
			drugDecay = 0.01;
		}

		for (int i = 0; i < length; i++){
			double randomValue = rn.Double();

			if (randomValue < ratioHealthy){
				Cells c = NewAgentSQ(i);
				c.CellInit(true,false, false, false);
			}
			else if(randomValue > ratioHealthy && randomValue < ratioHealthy + ratioInfected ) {
				Cells c = NewAgentSQ(i);
				c.CellInit(false,true, false, false);
			}
			else {
				Cells c = NewAgentSQ(i);
				c.CellInit(false,false, false, true);
			}
		}
	}


	public void ModelStep(int tick){

		for (Cells cell : this){
			if (cell.CellType == 1){ // infected cell
				double addedVirusConcentration = VirusSource(tick, cell);
				double currentVirusConcentration = virusConcentration.Get(cell.Isq());
				double newVirusConcentration = addedVirusConcentration + currentVirusConcentration;
				virusConcentration.Set(cell.Isq(), newVirusConcentration);
			}
		}

		virusConcentration.DiffusionADI(diffCoeff);
		virusConcentration.Update();

		double addedDrugConcentration = DrugSource(tick);
		for (Cells cell : this){
			//if (cell.CellType == 3){ // capillary i.e. drug source
				double currentDrugConcentration = drugConcentration.Get(cell.Isq());
				double newDrugConcentration = addedDrugConcentration + currentDrugConcentration;
				drugConcentration.Set(cell.Isq(), newDrugConcentration);

			//}
		}

		for (int i = 0; i < MAX_PDE_STEP; i++) {
			drugConcentration.DiffusionADI(diffCoeff);
			if (drugConcentration.MaxDelta() < threshold) {
				drugConcentration.Update();
				break;
			}

			drugConcentration.Update();
		}

		// virusConcentration decreases at each time step for all cells
		for (Cells cell : this){
			double virusNow = virusConcentration.Get(cell.Isq());
			double drugNow = drugConcentration.Get(cell.Isq());
			// double removalEfficacy = 2/(1+Math.exp(100*drugNow));
			// double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
			double drugVirusRemovalEff = 0; // nirmatrelvir has 0 effect on virus removal
			double totalVirusRemoval = virusRemovalRate + drugVirusRemovalEff;
			virusConcentration.Add(cell.Isq(), -totalVirusRemoval * virusNow);
		}
		virusConcentration.Update();

		// decay of the drug
		for (Cells cell : this){
			double drugNow = drugConcentration.Get(cell.Isq());
			drugConcentration.Add(cell.Isq(), - drugDecay * drugNow);
		}
		drugConcentration.Update();


		for (Cells cell : this){
			cell.CellInfection();
		}
		for (Cells cell : this) {
			cell.CellDeath();
		}
	}

	double TotalVirusConcentration() {
		double totalVirusConcentration = 0;
		for (int i = 0; i < length; i++){
			cellularVirusConcentration[i] = virusConcentration.Get(i);
		}
		for (double virusConcentrationInCell : cellularVirusConcentration ){
			totalVirusConcentration = totalVirusConcentration + virusConcentrationInCell;
		}
		return totalVirusConcentration;
	}


	double TotalDrugConcentration() {
		double totalDrugConcentration = 0;
		for (int i = 0; i < length; i++){
			cellularDrugConcentration[i] = drugConcentration.Get(i);
		}
		for (double drugConcentrationInCell : cellularDrugConcentration ){
			totalDrugConcentration = totalDrugConcentration + drugConcentrationInCell;
		}
		return totalDrugConcentration;
	}


	double VirusSource(int tick, Cells cell){

		double drugNow = drugConcentration.Get(cell.Isq());
		double drugVirusProdEff = 100 * Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
		return virusMax * (1-drugVirusProdEff);

	}

	double DrugSource(int tick){
		if (tick % 200 == 50) {
			return drugMax;
		} else {
			return 0.0;
		}
	}

	public void DrawModel(GridWindow vis){
		for (int i = 0; i < length; i++) {
			Cells drawMe = GetAgent(i);

			if (drawMe == null){
				vis.SetPix(i, RGB256(255, 255, 255));
			}
			else{
				if (drawMe.CellType == 0){		// Healthy cells
					vis.SetPix(i, RGB256(119, 198, 110));
				}
				else if (drawMe.CellType == 1){   // Infected cells
					vis.SetPix(i, RGB256(124, 65, 120));
				}
				else if (drawMe.CellType == 2){
					vis.SetPix(i, RGB(0, 0, 0));
				}
				else if (drawMe.CellType == 3){
					vis.SetPix(i, RGB(255, 255, 255));
				}
			}

			vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapRGB(virusConcentration.Get(i)));

			vis.SetPix(ItoX(i) + xDim*2, ItoY(i), HeatMapGBR(drugConcentration.Get(i)));
		}
	}

}


class Cells extends AgentSQ2Dunstackable<NirmatrelvirExperiments>{
	int CellType;

	public void CellInit(boolean isHealthy, boolean isInfected,
						 boolean isDead, boolean isCapillary){
		if(isHealthy == true){
			this.CellType = 0;
		}
		else if(isInfected == true){
			this.CellType = 1;
		}
		else if(isDead == true) {
			this.CellType = 2;
		}
		else if(isCapillary == true) {
			this.CellType = 3;
		}

	}

	public void CellInfection(){

		double drugConcentrationAtCell = G.drugConcentration.Get(Isq());
		double virusConAtCell = G.virusConcentration.Get(Isq());

		// double drugeff = 2/(1+Math.exp(100*chemoconcentr));
		// we consider a sigmoid function for drug efficacy
		// double drugInfectionRedEff = 100*Math.pow(drugConcentrationAtCell, 2)/(1+100*Math.pow(drugConcentrationAtCell,2));
		double drugInfectionRedEff = 0; // nirmatrelvir has no effect on infection rate
		double infecRate = G.infectionRate * (1 - drugInfectionRedEff) * G.xDim * G.yDim * virusConAtCell;

		if (this.CellType == 0){ // healthy cell
			if (G.rn.Double() < infecRate) {
				this.CellType = 1;
			}
		}
	}

	public void CellDeath(){

		if (this.CellType == 1) { // infected
			if(G.rn.Double() < G.deathProb){
				this.CellType = 2;
			}
		}
	}

}