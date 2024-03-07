package benders;

/*******
 * This section is the code of Benders Master Problem.
********/
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import ilog.cplex.IloCplex;
import common.Configure;
import common.Read_In;
import java.util.Arrays;
import ilog.concert.*;

public class BMP {
	public static String readPath = "data\\";
	public static int ord = 13;										//unused;
	public static ArrayList<ArrayList<Integer>> ARR = new ArrayList<ArrayList<Integer>>();
	public static int [] Machine_list_s = {1,2,3,5};				//unused;
	public static int [] Machine_list_m = {1,2,3,5};				//unused;
	public int m1;
	public int m2;
	public static double M = 999.0;
	public double run_time;
	public double makespan;
	public ArrayList<ArrayList<Integer>> solution;
	public String instance_name;
	public IloCplex cplex;
	public IloObjective obj;
	public double mpobj;
	
	//parameters;
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
	
	public BMP(String arr) {
		this.instance_name = arr;
	}
	
	
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
	public void model_constract(Read_In ri )throws IOException{
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
			
			//Annotate or add the following inequality to test its effectiveness
			//---constraints(6-7)---(Valid inequalities to strengthen the BMP)//
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
						int ESm = 0;
						for(int ii : this.N_0) {
							for(int kk : this.K_0) {
								if(this.ts[i][k][ii][kk] >  ESm)
									ESm = this.ts[i][k][ii][kk];
							}
						}
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr2 = this.cplex.sum(subexpr2,  ESm);
								}
							}
						}
					}
				}
				IloNumExpr subexpr3 = this.cplex.numExpr();
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addGe(subexpr3, this.C);
				cplex.addLe(subexpr1, this.C);
			}
			//---constraints(8-9)---(Valid inequalities to strengthen the BMP)//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int k : this.K) {
					subexpr1 = this.cplex.sum(subexpr1, this.Y[k][p]);
				}
				for(int k : this.K) {
					int ESm = 0;
					for(int kk : this.K_0) {
						if(this.tm[k][kk] >  ESm)
							ESm = this.tm[k][kk];
						}
					
					
					
					for(int h : this.K_0) {
						if(h != k) {
							subexpr2 = this.cplex.sum(subexpr2, 
									ESm);
						}
					}
				}
				IloNumExpr subexpr3 = this.cplex.numExpr();
				subexpr3 = this.cplex.sum(subexpr1, subexpr2);
				cplex.addGe(subexpr3, this.C);
				cplex.addLe(subexpr1, this.C);
			}
			
			

			
		}catch(IloException e) {
		System.err.println("Concert exception caught '" + e + "' caught");
		}
	}
	//solve the mode;
	public void model_solver()throws IOException{
		try {
				cplex.solve();
				if(!cplex.solve()) {
					System.exit(-1);
				}
				
				//output a bound of the problem;
				solution = new ArrayList<ArrayList<Integer>>();
				for(int p : this.Ps) {
					ArrayList<Integer> temp_sol = new ArrayList<Integer>();
					ArrayList<Integer> temp_n = new ArrayList<Integer>();
					ArrayList<Integer> temp_k = new ArrayList<Integer>();
					for(int i = 0; i < this.N_0.size(); i ++) {
						for(int j = 0; j <  this.K_0.size(); j ++) {
							if(Math.abs(cplex.getValue(this.XS[i][j][p]) - 1) < 0.00001) {
								temp_n.add(i); temp_k.add(j);
							}
						}
					}
					int sn = 0; int sk = 0;
					while(temp_n.size() > 0) {
						int remov = 0;
						for(int i = 0; i < temp_n.size(); i ++) {
							if(Math.abs(cplex.getValue(this.Lamd1[sn][sk][temp_n.get(i)][temp_k.get(i)][p]) - 1) < 0.00001) {
								remov = i;
								int job = 0;
								for(int k = 0; k < ARR.size(); k ++) {
									if(ARR.get(k).get(3) == temp_n.get(i)-1 && ARR.get(k).get(2) == temp_k.get(i)-1) {
										job = ARR.get(k).get(0);
									}
								}
								for(int j = 0; j < (cplex.getValue(this.Delta1[temp_n.get(i)][temp_k.get(i)][p]) ); j ++) {
									temp_sol.add(job);
								}
								sn = temp_n.get(i); sk = temp_k.get(i);
								break;
							}
						}
						temp_n.remove(remov); temp_k.remove(remov);
					}
					solution.add(temp_sol);
				}

				for(int p : this.Pm) {
					ArrayList<Integer> temp_sol = new ArrayList<Integer>();
					ArrayList<Integer> temp_k = new ArrayList<Integer>();
					for(int j = 0; j <  this.K_0.size(); j ++) {
						if(Math.abs(cplex.getValue(this.XM[j][p]) - 1) < 0.00001) {
							temp_k.add(j); 
						}
					}
					int sk = 0;
					while(temp_k.size() > 0) {
						int remov = 0;
						for(int i = 0; i < temp_k.size(); i ++) {
							if(Math.abs(cplex.getValue(this.Lamd2[sk][temp_k.get(i)][p]) - 1) < 0.00001) {
								remov = i;
								int job = 0;
								for(int n = 0; n <this.N_0.size(); n ++) {
									if(cplex.getValue(this.Delta2[n][temp_k.get(i)][p]) > 0.0001) {
										for(int k = 0; k < ARR.size(); k ++) {
											if(ARR.get(k).get(3) == n - 1 && ARR.get(k).get(2) == temp_k.get(i) - 1) {
												job = ARR.get(k).get(0);
											}
										}
										for(int j = 0; j < (cplex.getValue(this.Delta2[n][temp_k.get(i)][p]) ); j ++) {
											temp_sol.add(job);
										}
									}
								}
								sk = temp_k.get(i);
								break;
							}
						}
						 temp_k.remove(remov);
					}
					solution.add(temp_sol);
				}
				//output
				for(int i = 0; i < solution.size(); i ++) {
					for(int j = 0; j < solution.get(i).size(); j ++) {
						System.out.print(solution.get(i).get(j) + "  ");
					}
					System.out.println();
				}
				this.mpobj = cplex.getObjValue();
				System.out.println(cplex.getObjValue());
				
				
				
				
				File solution_file = new File("Solution_" + instance_name +".txt");
				FileWriter sol_writer = new FileWriter(solution_file);
				for(int ii = 0; ii < solution.size(); ii ++) {
					sol_writer.write("Machine: ");
					for(int jj = 0; jj <  solution.get(ii).size(); jj ++) {
						sol_writer.write( solution.get(ii).get(jj) + "  ");
					}
					sol_writer.write('\n');
				}
				double aaaa = cplex.getObjValue();
				String bbbb = String.valueOf(aaaa);
				sol_writer.write(bbbb);
				sol_writer.close();
				
				
				
				
				
				cplex.writeSolution("MP_relaxation.txt");
				makespan = cplex.getObjValue();
			}catch(IloException e) {
				System.err.println("Concert exception caught '" + e + "' caught");
		}
	}
}
