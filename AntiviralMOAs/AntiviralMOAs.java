package AntiviralMOAs;

// Experiment nr 1: reduced cell infection
// Experiment nr 2: reduced virus production (like nirmatrelvir)

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;
import HAL.Rand;
import static HAL.Util.*;

import java.io.File;

public class AntiviralMOAs{

    public static final int BIG_VALUE = (Integer.MAX_VALUE-1) / 2;

    public static void main(String[] args) {

        int y = 200, x = 200, visScale = 2;

        boolean isBoostedSlowerDecay = false;
        AppliedMOAsInExperiment appliedMOAsInExperiment = new AppliedMOAsInExperiment(0.0, 1800.0, 0.0, 0.0, 0.0, isBoostedSlowerDecay);

        GridWindow win = new GridWindow("Cellular state space, virus concentration.", x*2, y, visScale,true);

        NewExperiment experiment = new NewExperiment(x, y, visScale, new Rand(1), appliedMOAsInExperiment, 12*60, 0.2, 110.0);
        experiment.numberOfTicks = experiment.numberOfTicksDelay + experiment.numberOfTicksDrug;

        experiment.Init();
        double remainingHealthyCells = experiment.RunExperiment(win);


        System.out.println("In vivo MOA1 (infection) drug source [ng / ml]: " + experiment.drugs[0].drugSourceStomach);
        System.out.println("In vivo MOA2 (production) drug source [ng / ml]: " + experiment.drugs[1].drugSourceStomach);

        System.out.println("Remaining healthy cells: " + remainingHealthyCells);

        win.Close();

    }

}

class AppliedMOAsInExperiment {

    double drugSourceStomachReducedInfection = 0.0;
    double drugSourceStomachReducedProduction = 1800.0;
    double drugSourceStomachAcceleratedDeathRate = 0.0;
    double drugSourceStomachIncreasedClearance = 0.0;
    double drugSourceStomachMolecularTrap = 0.0;
    boolean isBoostedProductionBlocker = false; // like ritonavir

    public AppliedMOAsInExperiment(double drugSourceStomachReducedInfection, double drugSourceStomachReducedProduction, double drugSourceStomachAcceleratedDeathRate, double drugSourceStomachIncreasedClearance, double drugSourceStomachMolecularTrap, boolean isBoostedProductionBlocker){
        this.drugSourceStomachReducedInfection = drugSourceStomachReducedInfection;
        this.drugSourceStomachReducedProduction = drugSourceStomachReducedProduction;
        this.drugSourceStomachAcceleratedDeathRate = drugSourceStomachAcceleratedDeathRate;
        this.drugSourceStomachIncreasedClearance = drugSourceStomachIncreasedClearance;
        this.drugSourceStomachMolecularTrap = drugSourceStomachMolecularTrap;
        this.isBoostedProductionBlocker = isBoostedProductionBlocker;
    }
}

class Drug {

    // general properties
    double EC50 = 62; // in nM = nanoMolars, [nM] = 10^-9 [mol/L]; https://www.fda.gov/media/155050/download
    double molarMassDrug = 499.535;

    // in vivo properties
    double drugDecay = 0.013; // see the paxlovid_config.nb Mathematica notebook
    double drugSourceStomach = 0.0; // see the paxlovid_config.nb Mathematica notebook
    double drugDecayStomach = 0.015; // see the paxlovid_config.nb Mathematica notebook

    public Drug(){

    }

    public double DrugVirusProdEff(double drugNow){

        // double drugVirusProdEff = 7000 * Math.pow(drugNow, 2)/(1+7000*Math.pow(drugNow,2));
        // drugNow is the drug concentration in nanograms / ml
        // drugNow needs to be converted to [nM]s, as IC50 is given in [nM]s
        double drugNowInNanoMolars = NgPerMlToNanomolars(drugNow);
        double drugVirusProdEff = 1 / ( 1 + (EC50 / StochasticDrug(drugNowInNanoMolars)));

        return drugVirusProdEff;

    }

    double StochasticDrug(double drug){

        double stdDevOfGaussian = drug / 100;

        Rand random = new Rand();
        double stochasticDrug = random.Gaussian(drug, stdDevOfGaussian);
        stochasticDrug = stochasticDrug > 0.0 ? stochasticDrug : 0.0;

        return stochasticDrug;

    }

    double NgPerMlToNanomolars(double drugNow){
        return drugNow * Math.pow(10,3) / molarMassDrug;
    }

}

class DrugMOAReducedInfection extends Drug{
    // in vivo constructor
    public DrugMOAReducedInfection(AppliedMOAsInExperiment appliedMOAsInExperiment){
        this.drugSourceStomach = appliedMOAsInExperiment.drugSourceStomachReducedInfection;
    }
}

class DrugMOAReducedVirusProduction extends Drug{

    boolean isBoostedSlowerDecay;

    // in vivo constructor
    public DrugMOAReducedVirusProduction(AppliedMOAsInExperiment appliedMOAsInExperiment){

        this.drugSourceStomach = appliedMOAsInExperiment.drugSourceStomachReducedProduction;
        this.isBoostedSlowerDecay = appliedMOAsInExperiment.isBoostedProductionBlocker;

        if (isBoostedSlowerDecay == true){
            drugDecay = drugDecay / 3.0;
            drugDecayStomach = drugDecayStomach / 4.0;
            drugSourceStomach = drugSourceStomach * 3.8;
        }
    }
}

class DrugMOAAcceleratedDeathRate extends Drug{
    // in vivo constructor
    public DrugMOAAcceleratedDeathRate(AppliedMOAsInExperiment appliedMOAsInExperiment){
        this.drugSourceStomach = appliedMOAsInExperiment.drugSourceStomachAcceleratedDeathRate;
    }
}

class DrugMOAIncreasedVirusClearance extends Drug{
    // in vivo constructor
    public DrugMOAIncreasedVirusClearance(AppliedMOAsInExperiment appliedMOAsInExperiment){
        this.drugSourceStomach = appliedMOAsInExperiment.drugSourceStomachIncreasedClearance;
    }
}

class DrugMOAMolecularTrap extends Drug{
    // in vivo constructor
    public DrugMOAMolecularTrap(AppliedMOAsInExperiment appliedMOAsInExperiment){
        this.drugSourceStomach = appliedMOAsInExperiment.drugSourceStomachMolecularTrap;
    }
}

class NewExperiment extends AgentGrid2D<Cells>{

    public int x = 200;
    public int y = 200;
    public int visScale = 2;
    public int numberOfTicksDelay;
    public int numberOfTicksDrug;
    public int numberOfTicks;
    public PDEGrid2D virusCon;
    public PDEGrid2D immuneResponseLevel; // similar to interferon concentrations, but more generic
    public double[] drugCon = new double[5]; // 1: reducedInfection, 2: reduced Production, 3: accelerated death rate 4: immune stim, increased virus clearance 5: molecular trap
    public double[] drugConStomach = new double[5];
    public Rand rn;
    public double[] cellularVirusCon = new double[length];
    public double[] cellularImmuneResponseLevel = new double[length];

    public double fixedDamageRate;

    public double ratioHealthy = 0.9995, ratioInfected = 0.0005, ratioCapillaries = 0;

    // the following parameters are identical to those given in Table 2 and 3 for SARS-CoV-2 in the RSOS article
    public double virusRemovalRate = 1.67 * Math.pow(10,-3); // mu_V
    public double virusMax = 3.72 * Math.pow(10,-3); // f_{i,j}
    public double infectionRate = 1.01 * Math.pow(10,-7); // beta in the ODE
    public double deathProb = 7.02 * Math.pow(10,-4); // P_D
    public double virusDiffCoeff = 0.2; // D_V [sigma^2 / min]

    public Drug[] drugs = new Drug[5];

    public double immuneResponseDecay = 0.0005;
    public double immuneResponseDiffCoeff = 0.1;

    public double MAX_PDE_STEP = 1;
    public double threshold = 0.000001;

    public AppliedMOAsInExperiment appliedMOAsInExperiment;

    public FileIO outFile;
    public FileIO paramFile;
    public FileIO concentrationsFile;
    public String outputDir;

    public NewExperiment(int xDim, int yDim, int visScale, Rand rn, AppliedMOAsInExperiment appliedMOAsInExperiment, int numberOfTicksDelay, double virusDiffCoeff, double fixedDamageRate){

        super(xDim, yDim, Cells.class);
        this.x = xDim;
        this.y = yDim;
        this.visScale = visScale;
        this.numberOfTicksDelay = numberOfTicksDelay;
        this.virusDiffCoeff = virusDiffCoeff;
        this.rn = rn;

        this.appliedMOAsInExperiment = appliedMOAsInExperiment;

        this.fixedDamageRate = fixedDamageRate;

        this.drugs[0] = new DrugMOAReducedInfection(appliedMOAsInExperiment);
        this.drugs[1] = new DrugMOAReducedVirusProduction(appliedMOAsInExperiment);
        this.drugs[2] = new DrugMOAAcceleratedDeathRate(appliedMOAsInExperiment);
        this.drugs[3] = new DrugMOAIncreasedVirusClearance(appliedMOAsInExperiment);
        this.drugs[4] = new DrugMOAMolecularTrap(appliedMOAsInExperiment);

        this.numberOfTicksDrug = 5 * 24 * 60;

        virusCon = new PDEGrid2D(xDim, yDim);
        immuneResponseLevel = new PDEGrid2D(xDim, yDim);
        virusCon.Update();
        immuneResponseLevel.Update();

        outputDir = this.OutputDirectory();
        outFile = new FileIO(outputDir.concat("/").concat("Out").concat(".csv"), "w");
        paramFile = new FileIO(outputDir.concat("/").concat("Param").concat(".csv"), "w");
        concentrationsFile = new FileIO(outputDir.concat("/").concat("concentrations").concat(".csv"), "w");

    }

    public void Init(){

        WriteHeader();

        for (int i = 0; i < length; i++){
            double randomValue = rn.Double();

            if (randomValue <= ratioHealthy){
                Cells c = NewAgentSQ(i);
                c.CellInit(true,false, false, false);
            }
            else if(randomValue > ratioHealthy && randomValue <= ratioHealthy + ratioInfected ) {
                Cells c = NewAgentSQ(i);
                c.CellInit(false,true, false, false);
            }
            else {
                Cells c = NewAgentSQ(i);
                c.CellInit(false,false, false, true);
            }
        }
    }

    public double RunExperiment(GridWindow win){

        double[] cellCounts = CountCells();
        // System.out.println(cellCounts[0]+", " + cellCounts[1] + ", " + cellCounts[2]);

        for (int tick = 0; tick < this.numberOfTicks; tick ++){

            if ((numberOfTicksDelay == AntiviralMOAs.BIG_VALUE) && (fixedDamageRate <= 100) && (((cellCounts[1] + cellCounts[2]) / 40000) * 100 >= fixedDamageRate)){
                this.numberOfTicksDelay = tick - 1;
                this.numberOfTicks = this.numberOfTicksDelay + this.numberOfTicksDrug;
                System.out.println("Diffusion coeff.: " + this.virusDiffCoeff + ". Damage info: " + fixedDamageRate + " percent damage was found at tick " + numberOfTicksDelay + ".");
            }

            TimeStep(tick);
            DrawModel(win);

            if( tick > 0 && ( (tick % (24*60)) == 0 ))
                win.ToPNG(outputDir + "day" + Integer.toString(tick/(24*60)) + ".jpg");

            double totalVirusCon = TotalVirusCon();
            double totalImmuneResponseLevel = TotalImmuneResponseLevel();
            cellCounts = CountCells();
            concentrationsFile.Write(totalVirusCon + "," + totalImmuneResponseLevel + "," + drugCon[0] + "," + drugCon[1] + "," + drugCon[2] + "," + drugCon[3] + "," + drugCon[4] + "," + drugConStomach[0] + "," + drugConStomach[1] + "," + drugConStomach[2] + "," + drugConStomach[3] + "," + drugConStomach[4] + "\n");
            outFile.Write(tick +"," + cellCounts[0] + "," + cellCounts[1]+
                    "," + cellCounts[2] + "," + totalVirusCon + "," + drugCon + "," + drugConStomach + "\n");

        }

        CloseFiles();

        return cellCounts[0];

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

    void TimeStepImmune(int tick){

        // decay of the immuneResponseLevel
        for (Cells cell : this){
            double immuneResponseNow = immuneResponseLevel.Get(cell.Isq());
            immuneResponseLevel.Add(cell.Isq(), -immuneResponseDecay * immuneResponseNow);
        }
        immuneResponseLevel.Update();

        // immune response level increases
        for (Cells cell : this){
            if (cell.CellType == 1){ // infected cells produce interferon
                double addedImmuneResponseLevel = ImmuneResponseSource(tick, cell);
                double currentImmuneResponseLevel = immuneResponseLevel.Get(cell.Isq());
                double newImmuneResponseLevel = addedImmuneResponseLevel + currentImmuneResponseLevel;
                immuneResponseLevel.Set(cell.Isq(), newImmuneResponseLevel);
            }
        }

        immuneResponseLevel.DiffusionADI(immuneResponseDiffCoeff);
        immuneResponseLevel.Update();

    }

    void TimeStepVirus(){

        // decay of the virus
        for (Cells cell : this){
            // double removalEfficacy = 2/(1+Math.exp(100*drugNow));
            // double removalEfficacy = 100*Math.pow(drugNow, 2)/(1+100*Math.pow(drugNow,2));
            double drugVirusRemovalEff = 0.0 * drugCon[2]; // TODO: drugCon[2] virus removal
            double immuneVirusRemovalEff = 1 / (1 + 1/(Math.pow(immuneResponseLevel.Get(cell.Isq()),2)));
            // immuneVirusRemovalEff = 0.0;
            virusCon.Add(cell.Isq(), -virusRemovalRate * virusCon.Get(cell.Isq()));
            virusCon.Add(cell.Isq(), -drugVirusRemovalEff * virusCon.Get(cell.Isq()));
            virusCon.Add(cell.Isq(), -immuneVirusRemovalEff * virusCon.Get(cell.Isq()));
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

    void TimeStepDrug(int tick){

        for(int i = 0; i < 5; i++){

            // decay of the drug
            this.drugCon[i] -= this.drugs[i].drugDecay * this.drugCon[i];

            // decay of the drug in the stomach
            // and appearance of the drug at the lung epithelial cells
            double transferQuantity = this.drugs[i].drugDecayStomach * this.drugConStomach[i];
            this.drugCon[i] += transferQuantity;
            this.drugConStomach[i] -= transferQuantity;

            // drug appearance in the stomach
            this.drugConStomach[i] += DrugSourceStomach(i, tick);

        }

    }

    void TimeStepCells(){

        for (Cells cell : this){
            cell.CellInfection();
        }
        for (Cells cell : this) {
            cell.CellDeath();
        }

    }

    void TimeStep(int tick){

        TimeStepImmune(tick);
        TimeStepDrug(tick);
        TimeStepVirus();
        TimeStepCells();

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

    double TotalImmuneResponseLevel() {

        double totalImmuneResponseLevel = 0;
        for (int i = 0; i < length; i++){
            cellularImmuneResponseLevel[i] = immuneResponseLevel.Get(i);
        }
        for (double immuneResponseInCell : cellularImmuneResponseLevel ){
            totalImmuneResponseLevel = totalImmuneResponseLevel + immuneResponseInCell;
        }
        return totalImmuneResponseLevel;
    }

    double VirusSource(){

        return virusMax * (1 - this.drugs[1].DrugVirusProdEff(this.drugCon[1]));

    }

    double ImmuneResponseSource(int tick, Cells cell){

        return 0.0 * Math.pow(10,-3);

    }

    double DrugSourceStomach(int moa, int tick){

        if ((tick > numberOfTicksDelay)  && (((tick - numberOfTicksDelay) % (12 * 60)) == 1)) {
            return this.drugs[moa].drugSourceStomach;
        } else {
            return 0.0;
        }
    }

    void WriteHeader(){

        paramFile.Write("Parameters: \n init. ratio of healthy cells, virus removal rate, diff. coeff. \n");
        paramFile.Write(this.ratioHealthy + "," + this.virusRemovalRate + "," + this.virusDiffCoeff + "\n");
        outFile.Write("tick, healthy cells, infected cells, dead cells, "
                + "total virus conc., total drug conc., total drug conc. transfer \n");
    }

    void CloseFiles(){

        outFile.Close();
        paramFile.Close();
        concentrationsFile.Close();

    }

    String OutputDirectory(){

        java.util.Date now = new java.util.Date();
        java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
        String date_time = dateFormat.format(now);
        String projPath = PWD() + "/output/AntiviralMOAs";
        projPath += "/MOA1_" + appliedMOAsInExperiment.drugSourceStomachReducedInfection + "_MOA2_" + appliedMOAsInExperiment.drugSourceStomachReducedProduction;


        String outputDir = projPath + "/" + date_time + "__diff" + this.virusDiffCoeff;
        if(this.fixedDamageRate < 100.0){
            outputDir +=  "__damagerate" + this.fixedDamageRate + "/";
        } else {
            outputDir +=  "__delaytime" + this.numberOfTicksDelay + "/";
        }

        new File(outputDir).mkdirs();

        return outputDir;

    }

    void DrawModel(GridWindow vis){

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

            vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapRBG(virusCon.Get(i)));
        }
    }

}

class Cells extends AgentSQ2Dunstackable<NewExperiment>{
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

        double drugConAtCell = G.drugCon[0];
        double virusConAtCell = G.virusCon.Get(Isq());

        // we consider a sigmoid function for drug efficacy
        // double drugInfectionRedEff = 100*Math.pow(drugConAtCell, 2)/(1+100*Math.pow(drugConAtCell,2));
        double drugInfectionRedEff = drugConAtCell; // TODO
        double infectionProb = G.infectionRate * G.xDim * G.yDim; // converting beta to P_I (kind of)
        double effectiveInfectionProb = infectionProb * (1 - drugInfectionRedEff) * virusConAtCell;

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
