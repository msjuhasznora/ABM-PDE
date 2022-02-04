package covdrugInfecimmune;


//Spray drug and virus in a domain randomly. The drug affects on the 
//infection rate of the virus and also increase the virus removal rate

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Gui.GifMaker;
import HAL.Tools.FileIO;
import HAL.Rand;
//import LEARN_HERE.TumourClub.Cells;
import covdrugInfecimmune.Cells;
import covdrugInfecimmune.COVdrugbetaC;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import static HAL.Util.*;

public class COVdrugbetaC extends AgentGrid2D<Cells>{
public PDEGrid2D Virus;
public PDEGrid2D chemo;
public Rand rn;
public double[] Virusconcentr = new double[length];
public double[] chemoconcentr = new double[length];

public double Healthy = 0.998, Infected = 0.0002, 
Capillaries = 1-(Healthy+Infected) ;
public double consumrateVirus = 0.15, VirMax = 0.17, 
		chemoMax = 0.3;
public double infectionrate = 0.15 * Math.pow(10,-5);
public double deathprob = 0.011;
public double drugdecay = 0.004;
public double diffcoeff = 0.6; 
public double diffcoefch = 0.25;
//treatment parameters
public double perMed = 100;
public double drugperiod = 300;
double[] drugspace = new double[24];
double[][] chemoNow=new double[xDim][yDim];
public double MAX_PDE_STEP = 1;
public double Threshold = 0.000001;

public COVdrugbetaC(int xDim, int yDim, Rand rn){
   super(xDim, yDim, Cells.class);
   this.rn = rn;
   Virus = new PDEGrid2D(xDim, yDim);
   chemo = new PDEGrid2D(xDim, yDim);
   Virus.Update();
   chemo.Update();
}


public static void main(String[] args) {
   int y = 200, x = 200, visScale = 2;
   int Nt = 3200;//5000;
   int msPause = 0;


   //output file setup
   // TIME STAMP
   java.util.Date now = new java.util.Date();
   java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
   String date_time = dateFormat.format(now);
   // PATHS
   String projPath=PWD()+"/covdrugInfecimmune";
   String output_dir = projPath + "/output/" + date_time + "/";
   
   // CREATE OUTPUT DIR
	   new File(output_dir).mkdirs();
//   String output_dir = setting_dir;// + "/" + "MUTATION_PROB" + par_output ; \\Use with bash
//   new File(output_dir).mkdirs();


   //save data
   GridWindow win = new GridWindow("Model: cells, Virus, chemo", x*3, y, visScale,true);
   //FileIO outfile = new FileIO("setting_dir"+"Out" +".csv", "w");
   FileIO outfile = new FileIO(output_dir.concat("/").concat("Out").concat(".csv"), "w");
   FileIO paramfile = new FileIO(output_dir.concat("/").concat("Param").concat(".csv"), "w");

   COVdrugbetaC model =new COVdrugbetaC(x,y,new Rand(1));
   GifMaker gm=new GifMaker(output_dir.concat("/").concat("test").concat(".gif"),100,true);
//   win.ToGIF("test.jpg");

   //save parameter values to the file
	   paramfile.Write("Parameters: \n propHealthy, consumrateVirus, diffcoeff \n");
	   paramfile.Write(model.Healthy + "," + model.consumrateVirus + "," + model.diffcoeff + "\n");
	   outfile.Write("tick, HealthyCells, InfectedCells, DeadCells, "
		   		+ "BloodCells, VirusSummation, ChemoSummation \n");		
	   model.InitCov04(); // Put X For defining boundaries
	   model.DrawModel(win);
	   gm.AddFrame(win);

   double[] pops = model.CountCells();
   System.out.println(pops[0]+", " + pops[1] + ", " + pops[2]);

   for (int tick = 0; tick < Nt; tick ++){
       System.out.println(tick);
       model.ModelStep(tick);
       model.DrawModel(win);
//      gm.AddFrame(win);
       double summ = model.sum();
       double sumch = model.sumch();
	       double[] cellcount = model.CountCells();
	       outfile.Write(tick +"," + cellcount[0] + "," + cellcount[1]+
	    		   "," + cellcount[2]+"," + cellcount[3] + "," + summ+","+sumch+ "\n");


//       double[] cellcount = model.CountCells();
//       outfile.Write(tick +"," + cellcount[0] + "," + cellcount[1] + "," + cellcount[2] + "\n");

   }
   gm.Close();
   outfile.Close();
   paramfile.Close();
   win.Close();
}

double[] CountCells(){ //count each cell type
	   double HealthyCells = 0, InfectedCells = 0, DeadCells = 0,
			   CapillariesCells=0;
	   double[] cellcount = new double[4];
	   for (Cells cell: this){
	       if (cell.CellType == 0){
	       	HealthyCells+=1;
	       }
	       else if (cell.CellType == 1 ){
	       	InfectedCells+=1;
	       }
	       else if (cell.CellType == 2){
	       	DeadCells += 1;
	       }
	       else if (cell.CellType == 3){
	         	CapillariesCells += 1;
	         }

	    }
	    cellcount[0] = HealthyCells;
	    cellcount[1] = InfectedCells;
	    cellcount[2] = DeadCells;
	    cellcount[3] = CapillariesCells;


	   return cellcount;
	}

public void InitCov04(){
   //initiate the tumour: populate the grid with different cell types
	for (int i = 0; i < length; i++){
	       double randv = rn.Double();

	       if (randv < Healthy){
	           Cells c = NewAgentSQ(i);
	           c.CellInit(true,false, false, false);
	           //System.out.println("here  jglgh");
	       }
	       else if(randv > Healthy && randv < Healthy + Infected ) {
	           Cells c = NewAgentSQ(i);
	           c.CellInit(false,true, false, false);
	       }
	       else {
	      	 Cells c = NewAgentSQ(i);
	           c.CellInit(false,false, false, true);
	       }

	   }
   

	   //initial conditions for Virus 

// 	 double sigma = 15;
// 	 double meany = 100;
// 	 double meanx = 0;
//   double[][] res=new double[xDim][yDim];
//   double xx = 0;
//   double yy = 0;
//
//   for (int i = 0; i < xDim; i++) {
//   	for (int j = 0; j < yDim; j++) {
//   	
//   	res[i][j] =  (double) (Math.exp(-(Math.pow((xx-meanx),2)
//   			/(2*sigma*sigma))))*
//   			(Math.exp(-(Math.pow((yy-meany),2)/(2*sigma*sigma))));
//   	
//   	Virus.Set(i, j, res[i][j]);
//   	
//   	yy += 1;
//   	}
//   	yy = 0;
//   	xx = xx + 1;
//   }
   
   //initial condition for drug
   int jj = 0;
   for (int j = 0; j < yDim; j++){
   	for (int i = 0; i < xDim; i++) {
       chemo.Set(jj, chemoNow[j][i]);
       jj = jj + 1;
   	}
   }
   
//   for (int i = 0; i < 24; i++) {
//	   if(i == 0) {drugspace[i] = 0.68;}
//	   if(i == 1) {drugspace[i] = 0.29;}
//	   if(i == 2) {drugspace[i] = 0.75;}
//	   if(i == 3) {drugspace[i] = 0.43;}
//	   if(i == 4) {drugspace[i] = 0.20;}
//	   if(i == 5) {drugspace[i] = 0.55;}
//	   if(i == 6) {drugspace[i] = 0.36;}
//	   if(i == 7) {drugspace[i] = 0.08;}
//	   if(i == 8) {drugspace[i] = 0.36;}
//	   if(i == 9) {drugspace[i] = 0.88;}
//	   if(i == 10){drugspace[i] = 0.94;}
//	   if(i == 11){drugspace[i] = 0.26;}
//	   if(i == 12) {drugspace[i] = 0.82;}
//	   if(i == 13) {drugspace[i] = 0.45;}
//	   if(i == 14) {drugspace[i] = 0.35;}
//	   if(i == 15) {drugspace[i] = 0.09;}
//	   if(i == 16) {drugspace[i] = 0.42;} 
//	   if(i == 17) {drugspace[i] = 0.72;}
//	   if(i == 18) {drugspace[i] = 0.52;}
//	   if(i == 19) {drugspace[i] = 0.59;}
//	   if(i == 20) {drugspace[i] = 0.78;}
//	   if(i == 21) {drugspace[i] = 0.31;}
//	   if(i == 22) {drugspace[i] = 0.28;}
//	   if(i == 23) {drugspace[i] = 0.68;}
//	   if(i == 24) {drugspace[i] = 0.48;}
//   }
   

}
int ii = 0;
public void ModelStep(int tick){
   //Virus DYNAMICS

	double ViruConcentration = Diet(tick);
	   for (Cells cell : this){
		 if (cell.CellType == 1){ 
    	double Virnow = Virus.Get(cell.Isq());
    	double P = ViruConcentration + Virnow;
     	Virus.Set(cell.Isq(), P);
     
	  		 }
	   }

	   Virus.DiffusionADI(diffcoeff);
	   Virus.Update();

	   double DrugConcentration = ChemoCon(tick);
	   for (Cells cell : this){
		 if (cell.CellType == 3){ 
      	double Chemonow = chemo.Get(cell.Isq());
      	double Pd = DrugConcentration + Chemonow;
       	chemo.Set(cell.Isq(), Pd);
       
	  		 }
	   }

	   for (int i = 0; i < MAX_PDE_STEP; i++) {
		   chemo.DiffusionADI(diffcoeff);
//////	     double detla = Virus.MaxDelta();
	     if (chemo.MaxDelta() < Threshold) {
	    	 chemo.Update();
	    	 break;
	     }
	     
	     chemo.Update();
	 }

	   //Virus Decrease in each time step in all cells 
	   for (Cells cell : this){
	       double Virusnow = Virus.Get(cell.Isq());
	       double chemoconcentr = chemo.Get(cell.Isq());
//			double drugeff = 2/(1+Math.exp(100*chemoconcentr));
			double drugeff = 100*Math.pow(chemoconcentr, 2)/
					(1+100*Math.pow(chemoconcentr,2));
	       double immuneRes = consumrateVirus + drugeff;
	       Virus.Add(cell.Isq(), -immuneRes*Virusnow);
	   }
	   Virus.Update();
	
	      
	 
//if (tick == perMed ){
//	 
// double sigma = 30;
// 
//// double meanx = 200*rn.Double();
// //double meany = 200*rn.Double();
// 	 double meanx = 200*drugspace[ii];
//	 double meany = 200*drugspace[ii+1];
//	 ii += 1;
// perMed = perMed + drugperiod;
//    
// double xx = 0;
// double yy = 0;
//	  for (int i = 0; i < 200; i++) {
//		 for (int j = 0; j < 200; j++) {
//			   
//		  chemoNow[i][j] =  (double) (Math.exp(-(Math.pow((xx-meanx),2)/(2*sigma*sigma))))*
//			(Math.exp(-(Math.pow((yy-meany),2)/(2*sigma*sigma))));
//			   yy += 1;
//		   }
//		   yy = 0;
//		   xx = xx + 1;
//	   }
//	xx=0;
//int jj = 0;
//for (int i = 0; i < 200; i++){
//	for (int j = 0; j < 200; j++) {
////	double chemonow = chemo.Get(jj);
//	double updatechemo = chemoNow[i][j];
//   chemo.Add(jj,updatechemo);
//   jj = jj + 1;
//	}
//}
////chemo.Update();
//}
//   chemo.DiffusionADI(diffcoefch);
//   chemo.Update();

   
   //decay of the drug
   for (Cells cell : this){
	       double chemosnow = chemo.Get(cell.Isq());
	    chemo.Add(cell.Isq(), -drugdecay*chemosnow);
	   }
   chemo.Update();


	   //CELL STEP
//	   ShuffleAgents(rn);
	   for (Cells cell : this){
	//   	double VIRnow = Virus.Get(cell.Isq());
	   		cell.CellStep();
	   }
	   for (Cells cell : this) {
	   	cell.Celldeath();
	   }
}

double sum() {
	   double summ=0;
	   for (int i=0; i<length; i++){
	       Virusconcentr[i] = Virus.Get(i);
	   }
	   for (double num : Virusconcentr ){
	   	summ = summ + num;
	   }        	
	return summ;  
	}


double sumch() {
	  double summch=0;
	   for (int i=0; i<length; i++){
	       chemoconcentr[i] = chemo.Get(i);
	   }
	   for (double num : chemoconcentr ){
		 summch = summch + num;
	   }        	
	return summch;  
	}






double Diet(int tick){
   //Virus concentration in blood HealthyCells
   double VirusNow;

   //constant Virus
   VirusNow=VirMax;

//   if (tick < 5000){
//       VirusNow = VirMax;
//   }
//   else if (tick >= 100 && tick < 200){
//       VirusNow = VirMax;
//   }
//   else if (tick >= 200 && tick < 300){
//       VirusNow = VirMax;
//   }
//   else if (tick >= 300 && tick < 400){
//       VirusNow = VirMax;
//   }
//   else {
//       VirusNow = VirMax;
//   }

   return VirusNow;

}

double ChemoCon(int tick){
    //Virus concentration in blood HealthyCells
    double ChemoNow;

//    if (tick == perMed ){
//   	 perMed = perMed + drugperiod;
    
    //constant Virus
    ChemoNow=chemoMax;
//    }
//    else {
//   	 ChemoNow=0;
//    }
    
    
    
    return ChemoNow;

}

	public void DrawModel(GridWindow vis){
		   for (int i = 0; i < length; i++) {
		       Cells drawMe = GetAgent(i);
		
		       if (drawMe == null){
		           vis.SetPix(i, RGB256(255, 255, 255));
		       }
		       else{
		           if (drawMe.CellType == 0){		// Healthy cells are blue
		               vis.SetPix(i, RGB256(0, 255, 51));
		           }
		           else if (drawMe.CellType == 1){   // Infected cells are red
		               vis.SetPix(i, RGB256(255, 0, 0));
		           }
		           else if (drawMe.CellType == 2){
		               vis.SetPix(i, RGB(0, 0, 0));
		           }
		           else if (drawMe.CellType == 3){
		               vis.SetPix(i, RGB(255, 255, 255));
		           }
		       }
		
		       vis.SetPix(ItoX(i)+xDim, ItoY(i),HeatMapRGB(Virus.Get(i)));
		       
		       vis.SetPix(ItoX(i)+xDim*2, ItoY(i),HeatMapGBR(chemo.Get(i)));
		   }
		}

}




class Cells extends AgentSQ2Dunstackable<COVdrugbetaC>{
int CellType;

public void CellInit(boolean isVessel, boolean isInfectedCells,
		 boolean isDeadCells, boolean iscapillarCells){
   if(isVessel == true){
       this.CellType = 0;
   }
   else if(isInfectedCells== true){
       this.CellType = 1;
   }
   else if(isDeadCells== true) {
       this.CellType = 2;
   }
   else if(iscapillarCells== true) {
       this.CellType = 3;
   }

}

public void CellStep(){
		
		
		double chemoconcentr = G.chemo.Get(Isq());
//		double drugeff = 2/(1+Math.exp(100*chemoconcentr));
		double drugeff = 100*Math.pow(chemoconcentr, 2)/(1+100*Math.pow(chemoconcentr,2));
		// we consider sigmoid function for drug
		//double infectionrate = 0.3 * Math.pow(10,-5) * drugeff;
		double VIRnow = G.Virus.Get(Isq());
		double infecRate = G.infectionrate*(1-drugeff)*G.xDim*G.yDim* VIRnow;// No difference with the above VIRnow code
//		double p = G.infectionrate*G.xDim*G.yDim* VIRnow;
		if (this.CellType == 0){
			if (G.rn.Double() < infecRate) {
				this.CellType = 1;
			}
		}
	}

	public void Celldeath(){
	
	   
	   if (this.CellType == 1) {
	   	
	       if(G.rn.Double() < G.deathprob){
	       	this.CellType = 2;

	       }
	   }
	}


}
