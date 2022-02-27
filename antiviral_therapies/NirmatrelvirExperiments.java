package NirmatrelvirExperiments;

// based on Sadegh's CovDrugBetaC.java

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

	public PDEGrid2D virusCon;
	public PDEGrid2D drugCon;
	public Rand rn;
	public double[] cellularVirusCon = new double[length];
	public double[] cellularDrugCon = new double[length];

	public double ratioHealthy = 0.9995, ratioInfected = 0.0005, ratioCapillaries = 0;

	// the following parameters are identical to those given in Table 3 for SARS-CoV-2 in the RSOS article
	public double virusRemovalRate = 1.67 * Math.pow(10,-3); // mu_V
	public double virusMax = 1.66 * Math.pow(10,-2); // f_{i,j}
	public double infectionProb = 3.76 * Math.pow(10,-3); // P_I
	public double deathProb = 2.66 * Math.pow(10,-3); // P_D
	public double virusDiffCoeff = 0.2; // D_V [sigma^2 / min]

	public double drugMax = 5 * Math.pow(10,-2);
	public double drugDiffCoeff = 0.5;
	public double drugDecay = 3 * Math.pow(10,-3);

	public double MAX_PDE_STEP = 1;
	public double threshold = 0.000001;

	public boolean isRitonavirBoosted = true;

	public NirmatrelvirExperiments(int xDim, int yDim, Rand rn, boolean isRitonavirBoosted){
		super(xDim, yDim, Cells.class);
		this.rn = rn;
		this.isRitonavirBoosted = isRitonavirBoosted;
		virusCon = new PDEGrid2D(xDim, yDim);
		drugCon = new PDEGrid2D(xDim, yDim);
		virusCon.Update();
		drugCon.Update();
	}


	public static void main(String[] args) {
		int y = 200, x = 200, visScale = 2;
		int numberOfTicks = 5 * 24 * 60; // we follow the course of infection for 5 days, i.e. 5*24*60 minutes
		boolean isRitonavirBoosted = true;

		java.util.Date now = new java.util.Date();
		java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		String date_time = dateFormat.format(now);
		String projPath = PWD() + "/output/NirmatrelvirExperiments";
		if (isRitonavirBoosted == true){
			projPath += "/ritoBoosted";
		} else {
			projPath += "/nirmatrelvirOnly";
		}
		String outputDir = projPath + "/" + date_time + "/";
		new File(outputDir).mkdirs();

		GridWindow win = new GridWindow("Cellular state space, virus concentration, drug concentration", x*3, y, visScale,true);
		FileIO outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
		FileIO paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");

		NirmatrelvirExperiments model = new NirmatrelvirExperiments(x, y, new Rand(1), isRitonavirBoosted);
		GifMaker gm = new GifMaker(outputDir.concat("/").concat("test").concat(".gif"),100,true);
		// win.ToGIF("test.jpg");

		paramFile.Write("Parameters: \n init. ratio of healthy cells, virus removal rate, diff. coeff. \n");
		paramFile.Write(model.ratioHealthy + "," + model.virusRemovalRate + "," + model.virusDiffCoeff + "\n");
		outFile.Write("tick, healthy cells, infected cells, dead cells, "
				+ "total virus conc., total drug conc. \n");
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
			double totalVirusCon = model.TotalVirusCon();
			double totalDrugCon = model.TotalDrugCon();
			cellCounts = model.CountCells();
			outFile.Write(tick +"," + cellCounts[0] + "," + cellCounts[1]+
					"," + cellCounts[2] + "," + totalVirusCon + "," + totalDrugCon + "\n");

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
			drugDecay = 1 * Math.pow(10,-3);
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
				double addedVirusCon = VirusSource(tick, cell);
				double currentVirusCon = virusCon.Get(cell.Isq());
				double newVirusCon = addedVirusCon + currentVirusCon;
				virusCon.Set(cell.Isq(), newVirusCon);
			}
		}

		virusCon.DiffusionADI(virusDiffCoeff);
		virusCon.Update();

		double addedDrugCon = DrugSource(tick);
		for (Cells cell : this){
			//if (cell.CellType == 3){ // capillary i.e. drug source
				double currentDrugCon = drugCon.Get(cell.Isq());
				double newDrugCon = addedDrugCon + currentDrugCon;
				drugCon.Set(cell.Isq(), newDrugCon);

			//}
		}

		for (int i = 0; i < MAX_PDE_STEP; i++) {
			drugCon.DiffusionADI(drugDiffCoeff);
			if (drugCon.MaxDelta() < threshold) {
				drugCon.Update();
				break;
			}

			drugCon.Update();
		}

		// virusCon decreases at each time step for all cells
		for (Cells cell : this){
			double virusNow = virusCon.Get(cell.Isq());
			double drugNow = drugCon.Get(cell.Isq());
			// double removalEfficacy = 2/(1+Math.exp(100*drugNow));
			// double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
			double drugVirusRemovalEff = 0; // nirmatrelvir has 0 effect on virus removal
			double totalVirusRemoval = virusRemovalRate + drugVirusRemovalEff;
			virusCon.Add(cell.Isq(), -totalVirusRemoval * virusNow);
		}
		virusCon.Update();

		// decay of the drug
		for (Cells cell : this){
			double drugNow = drugCon.Get(cell.Isq());
			drugCon.Add(cell.Isq(), - drugDecay * drugNow);
		}
		drugCon.Update();


		for (Cells cell : this){
			cell.CellInfection();
		}
		for (Cells cell : this) {
			cell.CellDeath();
		}
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


	double TotalDrugCon() {
		double totalDrugCon = 0;
		for (int i = 0; i < length; i++){
			cellularDrugCon[i] = drugCon.Get(i);
		}
		for (double drugConInCell : cellularDrugCon ){
			totalDrugCon = totalDrugCon + drugConInCell;
		}
		return totalDrugCon;
	}


	double VirusSource(int tick, Cells cell){

		double drugNow = drugCon.Get(cell.Isq());
		double drugVirusProdEff = 100 * Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
		return virusMax * (1-drugVirusProdEff);

	}

	double DrugSource(int tick){
		if (tick % (12 * 60) == 50) {
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

			vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapRGB(virusCon.Get(i)));

			vis.SetPix(ItoX(i) + xDim*2, ItoY(i), HeatMapGBR(drugCon.Get(i)));
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

		double drugConAtCell = G.drugCon.Get(Isq());
		double virusConAtCell = G.virusCon.Get(Isq());

		// double drugeff = 2/(1+Math.exp(100*chemoconcentr));
		// we consider a sigmoid function for drug efficacy
		// double drugInfectionRedEff = 100*Math.pow(drugConAtCell, 2)/(1+100*Math.pow(drugConAtCell,2));
		double drugInfectionRedEff = 0; // nirmatrelvir has no effect on infection rate
		double effectiveInfectionProb = G.infectionProb * (1 - drugInfectionRedEff) * virusConAtCell;

		if (this.CellType == 0){ // healthy cell
			if (G.rn.Double() < effectiveInfectionProb) {
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