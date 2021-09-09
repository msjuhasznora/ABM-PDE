package covIteration1902;
//In this model the number of infected cells increases according to the virus concentration 
//At first there are some infected cells that diffuse virus in domain
//There is a probability that infected cells die and if virus concentration in each cell
//is more than a random variable that cell gets infected. 
//In initial step there is no infected cell and virus comes from a desire boundary condition
//in each cell stochastically
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Gui.UILabel;
import HAL.Gui.UIWindow;
import HAL.Gui.UIGrid;
import HAL.Gui.GifMaker;
import HAL.Tools.FileIO;
import HAL.Rand;
import HAL.Util;
//import LEARN_HERE.TumourClub.Cells;
import covIteration1902.Cells;
import covIteration1902.COV2Iterate;
import java.awt.Color;
import java.awt.Graphics;
//import java.awt.image.RenderedImage;
import java.io.File;
//import java.util.ArrayList;
//import java.util.Arrays;
import java.util.Random;
//import ColorBar.ColorBar;
//import javax.imageio.ImageIO;

import static HAL.Util.*;

public class COV2Iterate extends AgentGrid2D<Cells>{     // Agents
public PDEGrid2D Virus;
//public PDEGrid2D chemo;
public Rand rn;

//public int[] divHood = MooreHood(false);
public double[] Virusconcentr = new double[length];
//initial proportion of cell types
public  double propHealthy = 0.9995;

public double MAX_PDE_STEP = 1;
public double Threshold = 0.000001;

//public  double consumrateVirus = 2 * Math.pow(10, -2),
//VirMax = 4 * Math.pow(10, -2);
//public  double consumrateVirus = 2 * Math.pow(10, -2),
//		 VirMax = 4 * Math.pow(10, -2); //day
//public  double consumrateVirus = 2.08 * Math.pow(10, -3),
//VirMax = 1.666 * Math.pow(10, -2);

// COVID parameters
public  double consumrateVirus = 1.67 * Math.pow(10, -3),
VirMax = 3.72* Math.pow(10, -3);
public  double diffcoeff = 6*Math.pow(10,-1);
public double deathprob = 7.2 * Math.pow(10,-4);
public  double infectionrate = 1.01 * Math.pow(10,-7);
//public  double diffcoeff = 19 * Math.pow(10,-2);
//
//
//public  double infectionrate = 9.4 * Math.pow(10,-8);
//
//public double deathprob = 2.66 * Math.pow(10,-3);

//public  double diffcoeff = 2 * Math.pow(10,-1);
//public  double infectionrate = 1 * Math.pow(10,-6);
//public  double consumrateVirus = 2 * Math.pow(10, -2),
//		 VirMax = 4* Math.pow(10, -2);
//public double deathprob = 8 * Math.pow(10,-3);
//public static double[] iteration=new double[2000];
public static int Nt = 7000;
public static int Nofiter = 2000;
public static double[][] iteration=new double[Nofiter][Nt];
public static double[][] Viteration=new double[Nofiter][Nt];
//public static int jj = 0;
//	public double[] colors = new double[20];




public COV2Iterate(int xDim, int yDim, Rand rn){   // Constructor
   super(xDim, yDim, Cells.class);    // The last argument is a class type it uses  
   this.rn = rn;
   Virus = new PDEGrid2D(xDim, yDim);
//   chemo = new PDEGrid2D(xDim, yDim);
   Virus.Update();
//   chemo.Update();
}


public static void main(String[] args) {
	java.util.Date now = new java.util.Date();
	   java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
	   String date_time = dateFormat.format(now);
	   
	   // PATHS
	   String projPath=PWD()+"/covIteration1902";
	   final String output_dir = projPath + "/output/" + date_time + "/";
	   
	   // CREATE OUTPUT DIR
	   new File(output_dir).mkdirs();
	   FileIO VirusCon = new FileIO(output_dir.concat("/").concat("ViralLoad").concat(".csv"), "w");
	   FileIO outfile = new FileIO(output_dir.concat("/").concat("Out").concat(".csv"), "w");
//	   for (int t = 0; t < iteration[0].length ; t++){
//			   outfile.Write(t +",");
//			   VirusCon.Write(t +",");
//		 }
//	   outfile.Write("\n");
//	   VirusCon.Write("\n");
	   
	   // For iteration
	 for (int ii = 0; ii < Nofiter; ii++){
		 System.out.println(ii);
		 int y = 200, x = 200;
//		 visScale = 2;
//		 int msPause = 0;

		 COV2Iterate model =new COV2Iterate(x,y,new Rand(1));
   
		 model.InitCOV2Iterate(x); // Put X For defining boundaries
//		 model.DrawModel(win);
//   gm.AddFrame(win);
//		 win.ToPNG("Initial.jpg");

		 double[] pops = model.CountCells();
   	 	 System.out.println(pops[0]+", " + pops[1]+", " + pops[2]);

   	 	 for (int tick = 0; tick < Nt; tick ++){
//	   	  System.out.println(tick);
//	   	 	if(win.IsClosed()) {
//	   	 			 break;
//	   	 		 }
//	       win.TickPause(msPause);
	       model.ModelStep(tick,x);
//	       model.DrawModel(win);

	       double summ = model.sum();
	       double[] cellcount = model.CountCells();

	       iteration[ii][tick] = cellcount[3];
//	       Viteration[ii][tick] = summ;
	       
	       if (summ < 0.001) {
	    	   for (int t = tick; t < iteration[0].length-1 ; t++){
	    			 iteration[ii][t+1] = iteration[ii][t];
//	    			 Viteration[ii][t+1] = Viteration[ii][t];
	    	   	}
	    	   tick = Nt; 
	       }
   	 	 }
   	 	
//	   double[] pops2 = model.CountCells();
//	   System.out.println(pops2[0]+", " + pops2[1]+", " + pops2[2]);
//   gm.Close();

//   win.Close();
	 }
	 // To save data

	 for (int ii = 0; ii < iteration.length; ii++) {
		 for (int t = 0; t < iteration[0].length ; t++){
		  
			 outfile.Write(iteration[ii][t] +",");
//			 VirusCon.Write(Viteration[ii][t] +",");
		 }
		 outfile.Write("\n");
//		 VirusCon.Write("\n");
	}
	 outfile.Close();
//	 VirusCon.Close();
	 
//	 win.Close();

}


double[] CountCells(){ //count each cell type
   double HealthyCells = 0, InfectedCells = 0, DeadCells = 0, AfterInfected = 0; 
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
    	   AfterInfected += 1;
       }

   }
   cellcount[0] = HealthyCells;
   cellcount[1] = InfectedCells;
   cellcount[2] = DeadCells;
   cellcount[3] = AfterInfected;

   return cellcount;
}


public Random rnd = new Random();
double randomGenerator() {
   return rnd.nextDouble();
}


public void InitCOV2Iterate(int x){
	 //int j = 0;
   //initiate the tumour: populate the grid with different cell types
	double randv = randomGenerator();
	
	double R = Math.round(randv * 10000.0);
    R = R/10000.0;
   for (int i = 0; i < length; i++){

      
       if( i == Math.round(R * xDim*yDim)) {
           Cells c = NewAgentSQ(i);
//           c.CellInit(true,false,false);
           c.CellInit(false,true,false,false);
       }
       else {
//  	 else {
      	 //if (i>=x*(x-1)) {
           Cells c = NewAgentSQ(i);
//           c.CellInit(false,true,false);
           c.CellInit(true,false,false,false);
      	 //}
       }
       
   }

}

public void ModelStep(int tick,int xx){

   //update concentration within the InfectedCells
   double ViruConcentration = Diet(tick);
   for (Cells cell : this){
  	 
       if (cell.CellType == 1){ 
      	double Virnow = Virus.Get(cell.Isq());
      	double P = ViruConcentration + Virnow;
       	Virus.Set(cell.Isq(), P);
       }
   }
   int j = 0;
   for (int i = 0; i < MAX_PDE_STEP; i++) {
//// 
  	 j = j + 1;
   Virus.DiffusionADI(diffcoeff);
//////   double detla = Virus.MaxDelta();
   if (Virus.MaxDelta() < Threshold) {
  	 Virus.Update();
  	 break;
   }
   
   Virus.Update();
}


//   Virus.DiffusionADI(diffcoeff);
//   Virus.Update();
   
   //Virus Decrease in each time step in all cells 
   for (Cells cell : this){
       double Virusnow = Virus.Get(cell.Isq());
       Virus.Add(cell.Isq(), -consumrateVirus*Virusnow);
   }
   Virus.Update();

   //CELL STEP
   for (Cells cell : this){
   		cell.CellStep();
   }
   for (Cells cell : this) {
   	cell.Celldeath();
   }
}

double sum() {
   double summ=0;
   for (int i=0; i<length; i++){
       Virusconcentr[i] =	Virus.Get(i);
   }
   for (double num : Virusconcentr ){
   	summ = summ + num;
   }        	
return summ;  
}

double Diet(int tick){
   //glucose concentration in blood vessels
   double ViruNow;

   //constant glucose
   ViruNow=VirMax;

//   if (tick < 5000){
//       glucoseNow = gluMax;
//   }
//   else if (tick >= 100 && tick < 200){
//       glucoseNow = gluMin;
//   }
//   else if (tick >= 200 && tick < 300){
//       glucoseNow = gluMax;
//   }
//   else if (tick >= 300 && tick < 400){
//       glucoseNow = gluMin;
//   }
//   else {
//       glucoseNow = gluMax;
//   }

   return ViruNow;

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
//	    		 vis.SetPix(i, RGB256(185, 185, 185));
	    	 }
	    	 else if (drawMe.CellType == 1){   // Infected cells are red
	    		 vis.SetPix(i, RGB256(255, 0, 0));
	    	 }
	    	 else if (drawMe.CellType == 2){
//	    		 vis.SetPix(i, RGB256(128, 128, 128));
	    		 vis.SetPix(i, RGB256(0, 0, 0));
	    	 }
	    	 else if (drawMe.CellType == 3){
//	    		 vis.SetPix(i, RGB256(128, 128, 128));
	    		 vis.SetPix(i, RGB256(180, 3, 255));
	    	 }
	    	 vis.SetPix(ItoX(i)+xDim, ItoY(i),HeatMapRGB(Virus.Get(i)));
	    	 
	    	 	    	 
	    	// vis.SetPix(ItoX(i)+xDim, ItoY(i),HeatMapRGB(Virus.Get(i)));
//	    	 vis.SetPix(ItoX(i), ItoY(i),HeatMapGBR(Virus.Get(i)));
//       vis.SetPix(ItoX(i)+xDim*2, ItoY(i),HeatMapRGB(chemo.Get(i)));
	     }
   }

   
}


}





class Cells extends AgentSQ2Dunstackable<COV2Iterate>{
int CellType;

//CELL TYPE: 0 - Healthy, 1 - Infected, 2 - Deadcells

public void CellInit(boolean isHealthy, boolean isInfected,
		boolean isdead, boolean noeffect){
   if(isHealthy == true){
       this.CellType = 0;
   }
   else if(isInfected== true){
       this.CellType = 1;
   }
   else if(isdead== true) {
       this.CellType = 2;
   }
   else if(noeffect== true) {
       this.CellType = 3;
   }
}

public void CellStep(){
	
	 //double infectionrate =  1 * Math.pow(10,-6);
//	double infectionrate = 1.01 * Math.pow(10,-4); // Chentong article 
//	 double infectionrate = 0.94 * Math.pow(10,-7);
	 
	double VIRnow = G.Virus.Get(Isq());   // No difference with the above VIRnow code
//	double changestate = G.randomGenerator();
//	double q = 1 - (Math.pow(1-infectionrate, VIRnow));
	double p = G.infectionrate*G.xDim*G.yDim* VIRnow;
	if (G.randomGenerator() < p) {
//	if (G.rn.Double() < p) {//Virus concentration per each min
		// if virus concentration per minute is more than a random variable
		if (this.CellType == 0){
//			this.CellType = 1;
			this.CellType = 3; //to consider how much infected cells produce from one
		}
	}
}
public void Celldeath(){

	 
//   double deathprob = 8* Math.pow(10,-3);
//	 double deathprob = 2.66* Math.pow(10,-3);// !/5 times of infection rate according to the chinese article
//	 double deathprob = 0.11;
   if (this.CellType == 1) {
   	
       if(G.randomGenerator() < G.deathprob){
//  	 if(G.rn.Double() < G.deathprob){
       	this.CellType = 2;
//           this.Dispose();
//           return;
       }
   
   }
}


}

