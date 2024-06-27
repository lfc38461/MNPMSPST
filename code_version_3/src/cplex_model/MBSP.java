package cplex_model;

/*******
 * This section is the code of IMBT.
********/

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import ilog.cplex.IloCplex;
import common.Read_In;
import java.util.Arrays;
import java.util.List;
import common.Configure;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import ilog.concert.*;
public class MBSP {
	public static String readPath =  "data\\";
	public static int ord = 13;			//unused
	public static ArrayList<ArrayList<Integer>> ARR = new ArrayList<ArrayList<Integer>>();
	public int m1;
	public int m2;
	public static double M = 9999.0;
	public double run_time;
	public double makespan;
	
	public IloCplex cplex;
	public IloObjective obj;
	public double mip_ub;
	public double mip_lb;
	public double mip_gap;
	
	public ArrayList<Integer> J;
	public ArrayList<Integer> N;         
	public ArrayList<Integer> K;
	public ArrayList<Integer> J_0;
	public ArrayList<Integer> N_0;         
	public ArrayList<Integer> K_0;

	public ArrayList<Integer> Ps;
	public ArrayList<Integer> Pm;

	public int [] ws;
	public int [] wm;
	public int [][] ts;
	public int [][] tm;
	public int [][] FN;
	public int [][] FK;
	
	
	
	//variables set;
	public IloNumVar C;
	public IloNumVar [][][] XS;
	public IloNumVar [][]   XM;
	public IloNumVar [][] CS;
	public IloNumVar [][]   CM;
	public IloNumVar []   CPS;
	public IloNumVar []   CPM;
	public IloNumVar [][][]   YM;
	public IloNumVar [][]  OM;
	public IloNumVar [][][] Rho;
	
	
	
	//initialization
	public void initial(Read_In ri) {
		//---J---//
		this.J = new ArrayList<Integer>();
		int count = 0;
		for(int i = 0; i < ARR.size(); i ++) {
			count += ARR.get(i).get(6);
		}
		for(int i = 1; i < count + 1; i ++) {
			this.J.add(i);
		}
		//---J_0---//
		this.J_0 = new ArrayList<Integer>();
		for(int i = 0; i < count + 1; i ++) {
			this.J_0.add(i);
		}
		//---N---//
		this.N = new ArrayList<Integer>();
		for(int i = 1; i < ri.N + 1; i ++) {
			this.N.add(i);
		}
		//---K---//
		this.K = new ArrayList<Integer>();
		for(int i = 1; i < ri.K + 1; i ++) {
			this.K.add(i);
		}
		//---N_0---//
		this.N_0 = new ArrayList<Integer>();
		for(int i = 0; i < ri.N + 1; i ++) {
			this.N_0.add(i);
		}
		//---K_0---//
		this.K_0 = new ArrayList<Integer>();
		for(int i = 0; i < ri.K + 1; i ++) {
			this.K_0.add(i);
		}

		//---ps---//
		this.Ps = new ArrayList<Integer>();
		for(int i = 0; i < this.m1; i ++) {
			this.Ps.add(i);
		}
		//---pm---//
		this.Pm = new ArrayList<Integer>();
		for(int i = 0; i < this.m2; i ++) {
			this.Pm.add(i);
		}
		//---ws---//
		this.ws = new int [count + 1];
		this.ws[0] = 0;
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for(int i = 0; i < ARR.size(); i ++)
		{
			for(int j = 0; j < ARR.get(i).get(6); j ++) {
				temp.add(ARR.get(i).get(1));
			}
		}
		for(int i = 1; i <this.ws.length; i ++) {
			ws[i] = temp.get(i-1);
		}
		//---wm---//
		this.wm = new int [count + 1];
		this.wm[0] = 0;
		for(int i = 1; i <this.wm.length; i ++) {
			wm[i] = temp.get(i-1);
		}

		//---FN---//
		this.FN = new int [count + 1][ri.N + 1];
		for(int i = 0; i < count + 1; i ++) {	
			Arrays.fill(this.FN[i], 0);
		}
		int flg_0 = 0;
		for(int i = 0; i < ARR.size(); i ++) {
			for(int j = 0; j < ARR.get(i).get(6); j ++) {
				flg_0 ++;
				int a = ARR.get(i).get(3) + 1;
				this.FN[flg_0][a] = 1;
			}
		}
		
		//---FK---//
		this.FK = new int [count + 1][ri.K + 1];
		for(int i = 0; i < count + 1; i ++) {
			Arrays.fill(this.FK[i], 0);
		}
		int flg_1 = 0;
		for(int i = 0; i < ARR.size(); i ++) {
			for(int j = 0; j < ARR.get(i).get(6); j ++) {
				flg_1 ++;
				int b = ARR.get(i).get(2) + 1;
				this.FK[flg_1][b] = 1;
			}
		}
		
		//---ts---//
		this.ts = new int [count + 1][count + 1];
		ArrayList<Integer> temp_s = new ArrayList<Integer>();
		ArrayList<Integer> temp_i = new ArrayList<Integer>();
		for(int i = 0; i < ARR.size(); i ++)
		{
			for(int j = 0; j < ARR.get(i).get(6); j ++) {
				temp_s.add(ARR.get(i).get(4));
				temp_i.add(ARR.get(i).get(5));
			}
		}
		for(int i = 0; i < count + 1; i ++) {
			for(int j = 0; j < count + 1; j ++) {
				if(i == 0 && j != 0) {
					this.ts[i][j] = temp_i.get(j - 1);
				}
				else if(i != 0 && j != 0 && i != j) {
					int same = 0;
					for(int n = 0; n < ri.N + 1; n ++) {
						if(this.FN[i][n] == 1 && this.FN[j][n] == 1) {
							for(int k = 0; k < ri.K + 1; k ++) {
								if(this.FK[i][k] == 1 && this.FK[j][k] == 1) {
									same = 1;
								}
							}
						}
					}
					if(same == 0) {
						this.ts[i][j] = temp_s.get(j-1);
					}
					
				}
				else {
					this.ts[i][j] = 0;
				}
				
			}
		}
		//---tm---//
		this.tm = new int [ri.K + 1][ri.K + 1];
		for(int i = 0; i < ri.K + 1; i ++) {
			Arrays.fill(tm[i], 0);
		}
		for(int i = 0; i < ri.K + 1; i ++) {
			
			if(i == 0) {
				for(int k = 0; k < ARR.size(); k ++) {
					tm[i][ARR.get(k).get(2) + 1] = Math.max(ARR.get(k).get(5), tm[i][ARR.get(k).get(2) + 1]);
				}
			}
			else {
				for(int k = 0; k < ARR.size(); k ++) {
					tm[i][ARR.get(k).get(2) + 1] = Math.max(ARR.get(k).get(4), tm[i][ARR.get(k).get(2) + 1]);
				}
			}
			tm[i][i] = 0;
		}
		
		
	}
	
	//Bulid the MBSP model
	public void model_solver(Read_In ri )throws IOException{
		try {
			this.cplex = new IloCplex();
			this.cplex.setParam(IloCplex.Param.Simplex.Tolerances.Optimality, 1e-9);
			this.cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-9);
			this.cplex.setParam(IloCplex.DoubleParam.TimeLimit, 7200);
			this.cplex.setOut(null);
			
			//---C---//
			this.C = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "C");
			//---XS---//
			this.XS = new IloNumVar [ri.num + 1][ri.num + 1][m1];
			for(int i = 0; i < ri.num + 1; i ++) {
				for(int j = 0; j < ri.num + 1; j ++) {
							for(int kk = 0; kk < m1; kk ++) {
							XS[i][j][kk] = cplex.numVar(0, 1, IloNumVarType.Int, 
									"XS["+i+"]["+j+"]["+kk+"]");
							}
				}
			}
			//---XM---//
			this.XM = new IloNumVar [ri.num + 1][m2];
			for(int i = 0; i < ri.num + 1; i ++) {
					for(int ii = 0; ii < m2; ii ++) {
						XM[i][ii] = cplex.numVar(0, 1, IloNumVarType.Int, "XM["+i+"]["+ii+"]");
					}
			}
			//---CS---//
			this.CS = new IloNumVar [ri.num + 1][m1];
			for(int i = 0; i < ri.num + 1; i ++) {
				for(int ii = 0; ii < m1; ii ++) {
					CS[i][ii] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CS["+i+"]["+ii+"]");
				}
			}
			//---CM---//
			this.CM = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					CM[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CM["+i+"]["+j+"]");
				}
			}
			//---CPS---//
			this.CPS = new IloNumVar[m1];
			for(int i = 0; i < m1; i ++) {
				CPS[i] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CPS["+i+"]");
			}
			//---CPM---//
			this.CPM = new IloNumVar[m2];
			for(int i = 0; i < m2; i ++) {
				CPM[i] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CPM["+i+"]");
			}
			//---YM---//
			this.YM = new IloNumVar[ri.K + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						YM[i][j][ii] = cplex.numVar(0, 1, IloNumVarType.Int, "YM["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			//---OM---//
			this.OM = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					OM[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "OM["+i+"]["+j+"]");
				}
			}
			//---Rho---//
			this.Rho = new IloNumVar[ri.K + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						Rho[i][j][ii] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, 
								"Rho["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			
			
			
			
			//---OBJ---////---constraints(B.1)---//
			obj = cplex.addMinimize(C);
			
			//---constraints(B.2)---//
			for(int j : this.J) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				IloNumExpr subexpr3 = this.cplex.numExpr();
				for(int i : this.J_0) {
					for(int p : this.Ps) {
						if(j != i) {
							subexpr1 = this.cplex.sum(subexpr1, this.XS[i][j][p]);
						}
					}
				}
				for(int p : this.Pm) {
					subexpr2 = this.cplex.sum(subexpr2, this.XM[j][p]);
				}
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addEq(subexpr3, 1);
			}
			
			//---constraints(B.3)---//
			for(int i : this.J) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				IloNumExpr subexpr3 = this.cplex.numExpr();
				for(int j : this.J_0) {
					for(int p : this.Ps) {
						if(i != j) {
							subexpr1 = this.cplex.sum(subexpr1, this.XS[i][j][p]);
						}
					}
				}
				for(int p : this.Pm) {
					subexpr2 = this.cplex.sum(subexpr2, this.XM[i][p]);
				}
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addEq(subexpr3, 1);
			}
			
			//---constraints(B.4)---//
			for(int j : this.J) {
				for(int p : this.Ps) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int i : this.J_0) {
						subexpr1 = this.cplex.sum(subexpr1, this.XS[i][j][p]);
					}
					for(int h : this.J_0) {
						subexpr2 = this.cplex.sum(subexpr2, this.XS[j][h][p]);
					}
					cplex.addEq(subexpr1, subexpr2);
				}
			}
			
			//---constraints(B.5)---//
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int j : this.J) {
					subexpr1 = this.cplex.sum(subexpr1, this.XS[0][j][p]);
				}
				for(int j : this.J) {
					subexpr2 = this.cplex.sum(subexpr2, this.XS[j][0][p]);
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
			
			//---constraints(B.6)---//
			for(int i : this.J_0) {
				for(int j : this.J) {
					for(int p : this.Ps) {
						if(i != j) {
							IloNumExpr subexpr1 = this.cplex.numExpr();
							double subexpr2;
							subexpr1 = this.cplex.sum(
										this.cplex.diff(this.CS[j][p], this.CS[i][p])
										, this.cplex.prod(this.M, 
										this.cplex.diff(1, this.XS[i][j][p])));
							subexpr2 = this.ws[j] + this.ts[i][j];
							cplex.addGe(subexpr1, subexpr2);
						}
					}
				}
			}
			
			//---constraints(B.7)---//
			for(int k : this.K_0) {
				for(int h : this.K) {
					for(int p : this.Pm) {
						if(k != h) {
							IloNumExpr subexpr1 = this.cplex.numExpr();
							IloNumExpr subexpr2 = this.cplex.numExpr();
							subexpr1 = this.cplex.sum( this.cplex.diff(this.CM[h][p], this.CM[k][p]),
										this.cplex.prod(this.M, this.cplex.diff(1, this.YM[k][h][p])));
							subexpr2 = this.cplex.sum(this.tm[k][h], this.OM[h][p]);
							cplex.addGe(subexpr1, subexpr2);
						}
					}
				}
			}
			
			//---constraints(B.8)---//
			for(int a : this.N) {
				for(int k : this.K) {
					for(int p : this.Pm) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						for(int i : this.J) {
							IloNumExpr subexpr2 = this.cplex.numExpr();
							subexpr2 = this.cplex.prod( this.FK[i][k],
									   this.cplex.prod( this.FN[i][a], 
									   this.cplex.prod( this.wm[i], this.XM[i][p])));
							subexpr1 = this.cplex.sum(subexpr1, subexpr2);
						}
						cplex.addGe(this.OM[k][p], subexpr1);
					}
				}
			}
			
			//---constraints(B.9)---//
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				for(int i : this.J_0) {
					for(int j : this.J) {
						if(i != j ) {
							IloNumExpr subexpr2 = this.cplex.numExpr();
							double time = this.ts[i][j] + this.ws[j];
							subexpr2 = this.cplex.prod(this.XS[i][j][p], 
									   time);
							subexpr1 = this.cplex.sum(subexpr1, subexpr2);
						}
					}
				}
				cplex.addEq(subexpr1, this.CPS[p]);
			}
			
			//---constraints(B.10)-(B.11)--//
			for(int k : this.K_0) {
				for(int h : this.K) {
					for(int p : this.Pm) {
						if(k != h ) {
							IloNumExpr subexpr1 = this.cplex.numExpr();
							IloNumExpr subexpr2 = this.cplex.numExpr();
							IloNumExpr subexpr3 = this.cplex.numExpr();
							subexpr1 = this.cplex.sum(this.OM[h][p], this.tm[k][h]);
							subexpr2 = this.cplex.sum(subexpr1, 
									   		this.cplex.prod(this.M, 
											   this.cplex.diff(1, this.YM[k][h][p])));
							subexpr3 = this.cplex.diff(subexpr1, 
							   		this.cplex.prod(this.M, 
									   this.cplex.diff(1, this.YM[k][h][p])));
							cplex.addLe(this.Rho[k][h][p], subexpr2);
							cplex.addGe(this.Rho[k][h][p], subexpr3);
						}
					}
				}
			}
			
			//---constraints(B.12)---//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				for(int k : this.K_0) {
					for(int h : this.K) {
						if(k != h) {
							subexpr1 = this.cplex.sum(subexpr1, this.Rho[k][h][p]);
						}
					}
				}
				cplex.addEq(this.CPM[p], subexpr1);
			}
			
			//---constraints(B.13)---//
			for(int p : this.Ps) {
				cplex.addGe(this.C, this.CPS[p]);
			}
			
			
			//---constraints(B.14)---//
			for(int p : this.Pm) {
				cplex.addGe(this.C, this.CPM[p]);
			}
			
			//---constraints(B.15)---//
			for(int h : this.K) {
				for(int p : this.Pm) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int k : this.K_0) {
						if(k != h) {
							subexpr1 = this.cplex.sum(subexpr1, this.YM[k][h][p]);
							subexpr2 = this.cplex.sum(subexpr2, this.YM[h][k][p]);
						}
					}
					cplex.addEq(subexpr1, subexpr2);
					cplex.addLe(subexpr1, 1);
				}
			}
			
			//---constraints(B.16)---//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int k : this.K) {
					subexpr1 = this.cplex.sum(subexpr1, this.YM[0][k][p]);
					subexpr2 = this.cplex.sum(subexpr2, this.YM[k][0][p]);
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
			
			//---constraints(B.17)---//
			for(int k : this.K) {
				for(int i : this.J) {
					for(int p : this.Pm) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						IloNumExpr subexpr2 = this.cplex.numExpr();
						IloNumExpr subexpr3 = this.cplex.numExpr();
						subexpr1 = this.cplex.prod(this.M, 
								   this.cplex.diff(1, 
								   this.cplex.prod(this.XM[i][p], this.FK[i][k])));
						for(int h : this.K_0) {
							if(h != k ) {
								subexpr3 = this.cplex.sum(subexpr3, this.YM[k][h][p]);
							}
						}
						subexpr2 = this.cplex.diff(1, subexpr3);
						cplex.addGe(subexpr1, subexpr2);
					}
				}
			}
			
			
			
			
			cplex.solve();
			cplex.writeSolution("MP_3.txt");
			makespan = cplex.getObjValue();
			this.mip_ub = makespan;
			this.mip_gap = cplex.getMIPRelativeGap();
			this.mip_lb = this.mip_ub - (this.mip_ub * this.mip_gap);
			System.out.println(this.mip_lb);
			System.out.println(this.mip_ub);
			System.out.println(this.mip_gap);

			
		}catch(IloException e) {
		System.err.println("Concert exception caught '" + e + "' caught");
		}
	}
	
	
	
	
	public static void main(String[] args) throws IOException {	
		Configure conf = new Configure();
		conf.read_configure("configure\\configure.txt");
		
		ArrayList<ArrayList <String>> record = new ArrayList<ArrayList<String>>();
		ArrayList<String> csv_head = new ArrayList<String>();
		csv_head.add("Type");
		csv_head.add("Pattern");
		csv_head.add("Job");
		csv_head.add("Single");
		csv_head.add("Multi");
		csv_head.add("LB");
		csv_head.add("UB");
		csv_head.add("Gap");
		csv_head.add("RunTime");
		record.add(csv_head);
		for(int i = 0; i < conf.instance.size(); i ++) {
				Read_In read_in = new Read_In();
				read_in.read(readPath + conf.instance.get(i) + ".dta", ord, ARR);  


				MBSP model = new MBSP();
				model.m1 = Integer.parseInt(conf.sin_mac);
				model.m2 = Integer.parseInt(conf.mul_mac);
				System.out.println("Instance : " + conf.instance.get(i) + " M1 = " + model.m1+ " M2 = " + model.m2);
				
				model.initial(read_in);
				long starttime = System.currentTimeMillis();
				model.model_solver(read_in);
				long endtime = System.currentTimeMillis();
				model.run_time = (endtime-starttime)/1000.00000;
				System.out.println((endtime-starttime)/1000.00000 + "s");
				ARR.clear();
				
				ArrayList<String>record1 = new ArrayList<String>();
				record1.add(String.valueOf(read_in.K));
				record1.add(String.valueOf(4));
				record1.add(String.valueOf(read_in.num));
				record1.add(String.valueOf(model.m1));
				record1.add(String.valueOf(model.m2));
				record1.add(String.valueOf(model.mip_lb));
				record1.add(String.valueOf(model.mip_ub));
				record1.add(String.valueOf(model.mip_gap));
				record1.add(String.valueOf(model.run_time) + "s");
				record.add(record1);
				
		}
			
		
        FileWriter fileWriter = null;
        CSVPrinter csvPrinter = null;
        CSVFormat csvFormat = CSVFormat.DEFAULT;
		try {
		        fileWriter = new FileWriter("results\\MBSP_"+ conf.sin_mac +"_"+conf.mul_mac +".csv");
		        csvPrinter = new CSVPrinter(fileWriter, csvFormat);
		        for(int i = 0; i < record.size(); i ++) {
		        	csvPrinter.printRecord(record.get(i));
		        }
		        
				
				
			    } catch (Exception e) {
			        e.printStackTrace();
			    } finally {
			        try {
			            fileWriter.flush();
			            fileWriter.close();
			            csvPrinter.close();
			        } catch (Exception e) {
			            e.printStackTrace();
			        }
			    }

		
		}
	
	
}
