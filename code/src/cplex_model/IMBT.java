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


public class IMBT {
	//parameters;
	public static String readPath =  "data\\";
	public static int ord = 13;		//unused
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
	
	public ArrayList<Integer> N;         
	public ArrayList<Integer> K;
	public ArrayList<Integer> N_0;         
	public ArrayList<Integer> K_0;
	public int[][] q;
	public int[][] V;
	public ArrayList<Integer> Ps;
	public ArrayList<Integer> Pm;
	public ArrayList<Integer> P;
	public int [][] ws;
	public int [][] wm;
	public int [][][][] ts;
	public int [][] tm;
	
	//variables set;
	public IloNumVar C;
	public IloNumVar [][][] CS;
	public IloNumVar [][]   CM;
	public IloNumVar [][]   Y;
	public IloNumVar [][][] XS;
	public IloNumVar [][]   XM;
	public IloNumVar [][][] Delta1;
	public IloNumVar [][][] Delta2;
	public IloNumVar [][][][][] Lamd1;
	public IloNumVar [][][] Lamd2;
	public IloNumVar [][][] mu1;
	public IloNumVar [][]   mu2;
	public IloNumVar [][][] Rho1;
	public IloNumVar [][]   Rho2;
	
	
	
	//initialization;
	public void initial(Read_In ri) {
		this.N = new ArrayList<Integer>();
		for(int i = 1; i < ri.N + 1; i ++) {
			this.N.add(i);
		}
		this.K = new ArrayList<Integer>();
		for(int i = 1; i < ri.K + 1; i ++) {
			this.K.add(i);
		}
		this.N_0 = new ArrayList<Integer>();
		for(int i = 0; i < ri.N + 1; i ++) {
			this.N_0.add(i);
		}
		this.K_0 = new ArrayList<Integer>();
		for(int i = 0; i < ri.K + 1; i ++) {
			this.K_0.add(i);
		}
		//--index 0 is a virtual driection or class;
		this.q = new int [ri.N + 1][ri.K + 1];
		for(int i = 0; i < ARR.size(); i ++) {
			int a = ARR.get(i).get(3) + 1;
			int b = ARR.get(i).get(2) + 1;
			this.q[a][b] = ARR.get(i).get(6);
		}
		//---V---//
		this.V = new int[ri.N + 1][ri.K + 1];
		for(int i = 0; i < this.V.length; i ++) {
			Arrays.fill(this.V[i], 0);
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
		this.ws = new int [ri.N + 1][ri.K + 1];
		for(int i = 0; i < ARR.size(); i ++) {
			int a = ARR.get(i).get(3) + 1;
			int b = ARR.get(i).get(2) + 1;
			this.ws[a][b] = ARR.get(i).get(1);
		}
		//---wm---//
		this.wm = new int [ri.N + 1][ri.K + 1];
		for(int i = 0; i < ARR.size(); i ++) {
			int a = ARR.get(i).get(3) + 1;
			int b = ARR.get(i).get(2) + 1;
			this.wm[a][b] = ARR.get(i).get(1);
		}
		//---ts---//
		this.ts = new int [ri.N + 1][ri.K + 1][ri.N + 1][ri.K + 1];
		for(int i = 0; i < ri.N + 1; i ++) {
			for(int j = 0; j < ri.K + 1; j ++) {
				for(int k = 0; k < ri.N + 1; k ++) {
					Arrays.fill(ts[i][j][k], 0);
				}
				for(int k = 0; k < ARR.size(); k ++) {
					if(i == 0 && j == 0) {
						ts[i][j][ARR.get(k).get(3)+1][ARR.get(k).get(2)+1] = ARR.get(k).get(5);
					}
					else if(i != 0 && j != 0) {
						ts[i][j][ARR.get(k).get(3) + 1][ARR.get(k).get(2) + 1] = ARR.get(k).get(4);
					}
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
	
	//define the arc set V;
	public int [][] change_v(int [][] v, int a, int b ){
		for(int i = 0; i < v.length; i ++) {
			Arrays.fill(v[i], 0);
		}
		for(int i = 1; i < v.length; i ++) {
			for(int j = 1; j < v[i].length; j ++) {
				if((i != a) || (j != b)) {
					v[i][j] = 1;
				}
			}
		}
		v[0][0] = 1;
		return v;
	}
	//build the BMP model;
	public void model_solver(Read_In ri )throws IOException{
		try {
			this.cplex = new IloCplex();
			this.cplex.setParam(IloCplex.Param.Simplex.Tolerances.Optimality, 1e-9);
			this.cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-9);
			this.cplex.setParam(IloCplex.DoubleParam.TimeLimit, 7200);
			this.cplex.setOut(null);
			
			//---C---//
			this.C = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "C");;
			//---CS---//
			this.CS = new IloNumVar [ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m1; ii ++) {
						CS[i][j][ii] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CS["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			//---CM---//
			this.CM = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					CM[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "CM["+i+"]["+j+"]");
				}
			}
			//---Y---//
			this.Y = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					Y[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "Y["+i+"]["+j+"]");
					cplex.addGe(Y[i][j], 0);
				}
			}
			//---XS---//
			this.XS = new IloNumVar [ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
							for(int kk = 0; kk < m1; kk ++) {
							XS[i][j][kk] = cplex.numVar(0, 1, IloNumVarType.Int, 
									"XS["+i+"]["+j+"]["+kk+"]");
							cplex.addGe(XS[i][j][kk], 0);
							}
				}
			}
			//---XM---//
			this.XM = new IloNumVar [ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
					for(int ii = 0; ii < m2; ii ++) {
						XM[i][ii] = cplex.numVar(0, 1, IloNumVarType.Int, "XM["+i+"]["+ii+"]");
						cplex.addGe(XM[i][ii], 0);
					}
			}
			//---Delta---//
			this.Delta1 = new IloNumVar[ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m1; ii ++) {
						Delta1[i][j][ii] = cplex.numVar(0, q[i][j], IloNumVarType.Int, "Delta1["+i+"]["+j+"]["+ii+"]");
						cplex.addGe(Delta1[i][j][ii], 0);
					}
				}
			}
			this.Delta2 = new IloNumVar[ri.N + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						Delta2[i][j][ii] = cplex.numVar(0, q[i][j], IloNumVarType.Int, "Delta2["+i+"]["+j+"]["+ii+"]");
						cplex.addGe(Delta2[i][j][ii], 0);
					}
				}
			}
			//---lamd1---//
			this.Lamd1 = new IloNumVar [ri.N + 1][ri.K + 1][ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < ri.N + 1; ii ++) {
						for(int jj = 0; jj < ri.K + 1; jj ++) {
							for(int kk = 0; kk < m1; kk ++) {
								Lamd1[i][j][ii][jj][kk] = cplex.numVar(0, 1, IloNumVarType.Int, 
									"Lamd1["+i+"]["+j+"]["+ii+"]["+jj+"]["+kk+"]");
								cplex.addGe(Lamd1[i][j][ii][jj][kk], 0);
							}
						}
					}
				}
			}
			//---Lamd2---//
			this.Lamd2 = new IloNumVar [ri.K + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						Lamd2[i][j][ii] = cplex.numVar(0, 1, IloNumVarType.Int, "Lamd2["+i+"]["+j+"]["+ii+"]");
						cplex.addGe(Lamd2[i][j][ii], 0);
					}
				}
			}
			
			
			//---OBJ---////---constraints(1a)---//
			obj = cplex.addMinimize(C);

			
			//---constraints(1b)---//
			for(int i : this.N) {
				for(int k : this.K) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int p : Ps) {
						subexpr1 = this.cplex.sum(subexpr1, this.Delta1[i][k][p]);
					}
					for(int p : Pm) {
						subexpr2 = this.cplex.sum(subexpr2, this.Delta2[i][k][p]);
					}
					IloNumExpr subexpr3 = this.cplex.numExpr();
					subexpr3 = this.cplex.sum(subexpr1, subexpr2);
					cplex.addEq(subexpr3, q[i][k]);
				}
			}
					
			//---constraints(1c)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						IloNumExpr subexpr2 = this.cplex.numExpr();
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr1 = this.cplex.sum(subexpr1, this.Lamd1[h][g][i][k][p]);
									subexpr2 = this.cplex.sum(subexpr2, this.Lamd1[i][k][h][g][p]);
								}
							}
						}
						cplex.addEq(subexpr1,subexpr2);
						cplex.addLe(subexpr1, 1);
						cplex.addLe(subexpr2, 1);
					}
				}
			}
			//---constraints(1d)---//
			for(int k : this.K) {
				for(int p : this.Pm) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int h : this.K_0) {
						if(h != k) {
							subexpr1 = this.cplex.sum(subexpr1,this.Lamd2[h][k][p]);
							subexpr2 = this.cplex.sum(subexpr2,this.Lamd2[k][h][p]);
						}
					}
					cplex.addEq(subexpr1,subexpr2);
					cplex.addLe(subexpr1, 1);
					cplex.addLe(subexpr2, 1);
				}
			}
			//---constraints(1e)---//
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int i : this.N) {
					for(int k : this.K) {
						subexpr1 = this.cplex.sum(subexpr1, this.Lamd1[0][0][i][k][p]);
						subexpr2 = this.cplex.sum(subexpr2, this.Lamd1[i][k][0][0][p]);
					}
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
			
			//---constraints(1f)---//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int k : this.K) {
						subexpr1 = this.cplex.sum(subexpr1, this.Lamd2[0][k][p]);
						subexpr2 = this.cplex.sum(subexpr2, this.Lamd2[k][0][p]);
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
		
			//---constraints(1g)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						IloNumExpr subexpr2 = this.cplex.numExpr();
						subexpr2 = this.cplex.prod(this.M, this.XS[i][k][p]);
						cplex.addLe(this.Delta1[i][k][p], subexpr2);
					}
				}
			}
			
			//---constraints(1h)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Pm) {
						IloNumExpr subexpr2 = this.cplex.numExpr();
						subexpr2 = this.cplex.prod(this.M, this.XM[k][p]);
						cplex.addLe(this.Delta2[i][k][p], subexpr2);
					}
				}
			}
			
			//---constraints(1i)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									IloNumExpr subexpr1 = this.cplex.numExpr();
									IloNumExpr subexpr2 = this.cplex.numExpr();
									IloNumExpr subexpr3 = this.cplex.numExpr();
									IloNumExpr subexpr4 = this.cplex.numExpr();
									subexpr1 = this.cplex.diff(this.CS[h][g][p], this.CS[i][k][p]);
									subexpr2 = this.cplex.prod(this.M, this.cplex.diff(1, this.Lamd1[i][k][h][g][p]));
									subexpr3 = this.cplex.sum(this.ts[i][k][h][g], this.cplex.prod(this.ws[h][g], 
											   this.Delta1[h][g][p]));
									subexpr4 = this.cplex.sum(subexpr1, subexpr2);
									cplex.addGe(subexpr4, subexpr3);
								}
							}
						}
					}
				}
			}
			
			//---constraints(1j)---//
			for(int k : this.K) {
				for(int p : this.Pm) {
					for(int h : this.K) {
							if(h != k) {
								IloNumExpr subexpr1 = this.cplex.numExpr();
								IloNumExpr subexpr2 = this.cplex.numExpr();
								IloNumExpr subexpr3 = this.cplex.numExpr();
								IloNumExpr subexpr4 = this.cplex.numExpr();
								subexpr1 = this.cplex.diff(this.CM[h][p], this.CM[k][p]);
								subexpr2 = this.cplex.prod(this.M, this.cplex.diff(1, this.Lamd2[k][h][p]));
								subexpr3 = this.cplex.sum(this.tm[k][h], this.Y[h][p]);
								subexpr4 = this.cplex.sum(subexpr1, subexpr2);
								cplex.addGe(subexpr4, subexpr3);
							}
						}
					}
				}

			//---constraints(1k)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Pm) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						subexpr1 = this.cplex.prod(this.wm[i][k], this.Delta2[i][k][p]);
						cplex.addLe(subexpr1, this.Y[k][p]);
					}
				}
			}
	
			//---constraints(1l-1o)
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						IloNumExpr subexpr2 = this.cplex.numExpr();
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr1 = this.cplex.sum(subexpr1, this.Lamd1[i][k][h][g][p]);
									subexpr2 = this.cplex.sum(subexpr2, this.Lamd1[h][g][i][k][p]);
								}
							}
						}
						cplex.addEq(this.XS[i][k][p], subexpr1);
						cplex.addEq(this.XS[i][k][p], subexpr2);
					}
				}
			}
			for(int k : this.K) {
				for(int p : this.Pm) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int h : this.K_0) {
						if(h != k) {
							subexpr1 = this.cplex.sum(subexpr1, this.Lamd2[k][h][p]);
							subexpr2 = this.cplex.sum(subexpr2, this.Lamd2[h][k][p]);
						}
					}
					cplex.addEq(this.XM[k][p], subexpr1);
					cplex.addEq(this.XM[k][p], subexpr2);
				}
			}					
			
			//---constraints(1p)---//
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int i : this.N) {
					for(int k : this.K) {
						subexpr1 = this.cplex.sum(subexpr1, this.cplex.prod(this.ws[i][k], this.Delta1[i][k][p]));
					}
				}
				for(int i : this.N_0) {
					for(int k : this.K_0) {
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr2 = this.cplex.sum(subexpr2, this.cplex.prod
											(this.Lamd1[i][k][h][g][p], this.ts[i][k][h][g]));
								}
							}
						}
					}
				}
				IloNumExpr subexpr3 = this.cplex.numExpr();
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addLe(subexpr3, this.C);
			}
			
			//---constraints(1q)---//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int k : this.K) {
					subexpr1 = this.cplex.sum(subexpr1, this.Y[k][p]);
				}
				for(int k : this.K) {
					for(int h : this.K_0) {
						if(h != k) {
							subexpr2 = this.cplex.sum(subexpr2, this.cplex.prod(this.Lamd2[k][h][p],
									this.tm[h][k]));
						}
					}
				}
				IloNumExpr subexpr3 = this.cplex.numExpr();
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addLe(subexpr3, this.C);
			}
			
			
			cplex.solve();
			cplex.writeSolution("MP_2.txt");
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
		

	
	//run only for IMBT model;
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


				IMBT model = new IMBT();
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
		        fileWriter = new FileWriter("results\\IMBT_"+ conf.sin_mac +"_"+conf.mul_mac +".csv");
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
