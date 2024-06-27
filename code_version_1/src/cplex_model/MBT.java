package cplex_model;

/*******
 * This section is the code of MBT.
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

public class MBT {
	//parameters;
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
	public IloNumVar Z;
	public IloNumVar [][][] ZS;
	public IloNumVar [][]   ZM;
	public IloNumVar [][]   Y;
	public IloNumVar [][][][][] XS;
	public IloNumVar [][][]     XM;
	public IloNumVar [][][]     Delta1;
	public IloNumVar [][][]     Delta2;
	public IloNumVar [][][][][] FS; 	//rho
	public IloNumVar [][][] FM;			//rho       
	
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
	//build the MBT model;
	public void model_solver(Read_In ri )throws IOException{
		try {
			this.cplex = new IloCplex();
			this.cplex.setParam(IloCplex.Param.Simplex.Tolerances.Optimality, 1e-9);
			this.cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-9);
			this.cplex.setParam(IloCplex.DoubleParam.TimeLimit, 7200);
			this.cplex.setOut(null);
			
			//---Z---//
			this.Z = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "Z");;
			//---ZS---//
			this.ZS = new IloNumVar [ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m1; ii ++) {
						ZS[i][j][ii] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "ZS["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			//---ZM---//
			this.ZM = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					ZM[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "ZM["+i+"]["+j+"]");
				}
			}
			//---Y---//
			this.Y = new IloNumVar[ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < m2; j ++) {
					Y[i][j] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "Y["+i+"]["+j+"]");
				}
			}
			//---XS---//
			this.XS = new IloNumVar [ri.N + 1][ri.K + 1][ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < ri.N + 1; ii ++) {
						for(int jj = 0; jj < ri.K + 1; jj ++) {
							for(int kk = 0; kk < m1; kk ++) {
							XS[i][j][ii][jj][kk] = cplex.numVar(0, 1, IloNumVarType.Int, 
									"XS["+i+"]["+j+"]["+ii+"]["+jj+"]["+kk+"]");
							}
						}
					}
				}
			}
			//---XM---//
			this.XM = new IloNumVar [ri.K + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						XM[i][j][ii] = cplex.numVar(0, 1, IloNumVarType.Int, "XM["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			//---Delta---//
			this.Delta1 = new IloNumVar[ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m1; ii ++) {
						Delta1[i][j][ii] = cplex.numVar(0, q[i][j], IloNumVarType.Int, "Delta1["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			this.Delta2 = new IloNumVar[ri.N + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						Delta2[i][j][ii] = cplex.numVar(0, q[i][j], IloNumVarType.Int, "Delta2["+i+"]["+j+"]["+ii+"]");
					}
				}
			}

			//---FS---//
			this.FS = new IloNumVar [ri.N + 1][ri.K + 1][ri.N + 1][ri.K + 1][m1];
			for(int i = 0; i < ri.N + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < ri.N + 1; ii ++) {
						for(int jj = 0; jj < ri.K + 1; jj ++) {
							for(int kk = 0; kk < m1; kk ++) {
							FS[i][j][ii][jj][kk] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Int, 
									"FS["+i+"]["+j+"]["+ii+"]["+jj+"]["+kk+"]");
							}
						}
					}
				}
			}
			//---FM---//
			this.FM = new IloNumVar [ri.K + 1][ri.K + 1][m2];
			for(int i = 0; i < ri.K + 1; i ++) {
				for(int j = 0; j < ri.K + 1; j ++) {
					for(int ii = 0; ii < m2; ii ++) {
						FM[i][j][ii] = cplex.numVar(0.0, Double.MAX_VALUE, IloNumVarType.Float, "FM["+i+"]["+j+"]["+ii+"]");
					}
				}
			}
			
			
			
			
			
			//---OBJ---////---constraints(1a)---//
			obj = cplex.addMinimize(Z);

			
			//---constraints(A.1b)---//
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
			//---constraints(A.1c)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						IloNumExpr subexpr2 = this.cplex.numExpr();
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr1 = this.cplex.sum(subexpr1, this.XS[h][g][i][k][p]);
									subexpr2 = this.cplex.sum(subexpr2, this.XS[i][k][h][g][p]);
								}
							}
						}
						cplex.addEq(subexpr1,subexpr2);
						cplex.addLe(subexpr1, 1);
						cplex.addLe(subexpr2, 1);
					}
				}
			}
			//---constraints(A.1d)---//
			for(int k : this.K) {
				for(int p : this.Pm) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					for(int h : this.K_0) {
						if(h != k) {
							subexpr1 = this.cplex.sum(subexpr1,this.XM[h][k][p]);
							subexpr2 = this.cplex.sum(subexpr2,this.XM[k][h][p]);
						}
					}
					cplex.addEq(subexpr1,subexpr2);
					cplex.addLe(subexpr1, 1);
					cplex.addLe(subexpr2, 1);
				}
			}
			//---constraints(A.1e)---//
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int i : this.N) {
					for(int k : this.K) {
						subexpr1 = this.cplex.sum(subexpr1, this.XS[0][0][i][k][p]);
						subexpr2 = this.cplex.sum(subexpr2, this.XS[i][k][0][0][p]);
					}
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
			
			//---constraints(A.1f)---//
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				IloNumExpr subexpr2 = this.cplex.numExpr();
				for(int k : this.K) {
						subexpr1 = this.cplex.sum(subexpr1, this.XM[0][k][p]);
						subexpr2 = this.cplex.sum(subexpr2, this.XM[k][0][p]);
				}
				cplex.addEq(subexpr1, 1);
				cplex.addEq(subexpr2, 1);
			}
			
			
			//---constraints(A.1g)---//
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
									IloNumExpr subexpr5 = this.cplex.numExpr();
									IloNumExpr subexpr6 = this.cplex.numExpr();
									subexpr1 = this.cplex.diff(this.ZS[i][k][p], this.ZS[h][g][p]);
									subexpr2 = this.cplex.diff(1, this.XS[h][g][i][k][p]);
									subexpr3 = this.cplex.prod(this.M, subexpr2);
									subexpr4 = this.cplex.sum(subexpr1, subexpr3);
									
									subexpr5 = this.cplex.prod(this.ws[i][k], this.Delta1[i][k][p]);
									subexpr6 = this.cplex.sum(this.ts[h][g][i][k], subexpr5);
									cplex.addGe(subexpr4, subexpr6);
									
								}
							}
						}
					}
				}
			}
			//---constraints(A.1h)---//
			for(int k : this.K) {
				for(int p : this.Pm) {
					for(int h : this.K_0) {
						if(h != k) {
							IloNumExpr subexpr1 = this.cplex.numExpr();
							IloNumExpr subexpr2 = this.cplex.numExpr();
							IloNumExpr subexpr3 = this.cplex.numExpr();
							subexpr1 = this.cplex.diff(this.ZM[k][p], this.ZM[h][p]);
							subexpr3 = this.cplex.diff(1, this.XM[h][k][p]);
							subexpr3 = this.cplex.prod(this.M, subexpr3);
							subexpr1 = this.cplex.sum(subexpr1, subexpr3);
							subexpr2 = this.cplex.sum(this.tm[h][k], this.Y[k][p]);
							cplex.addGe(subexpr1, subexpr2);
						}
					}
				}
			}
			
			
			//---constraints(A.1i)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Pm) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						subexpr1 = this.cplex.prod(this.wm[i][k], this.Delta2[i][k][p]);
						cplex.addGe(this.Y[k][p], subexpr1);
					}
				}
			}
			
			//---constraints(A.1j)---//
			IloNumExpr subexprl = this.cplex.numExpr();
			IloNumExpr subexprr = this.cplex.numExpr();
			for(int p : this.Ps) {
				subexprl = this.cplex.sum(subexprl,this.ZS[0][0][p]);
			}
			for(int p : this.Pm) {
				subexprr = this.cplex.sum(subexprr,this.ZM[0][p]);
			}
			cplex.addEq(subexprl, subexprr);
			cplex.addEq(subexprl, 0);
			cplex.addEq(subexprr, 0);
			
			//---constraints(1k)---//
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
									subexpr1 = this.cplex.prod(this.M, this.XS[h][g][i][k][p]);
									subexpr2 = this.cplex.sum(this.ts[h][g][i][k], this.cplex.prod(this.ws[i][k], this.Delta1[i][k][p]));
									subexpr3 = this.cplex.diff(subexpr2, this.cplex.prod(this.M, this.cplex.diff(1, this.XS[h][g][i][k][p])));
									cplex.addLe(FS[h][g][i][k][p], subexpr1);
									cplex.addLe(FS[h][g][i][k][p], subexpr2);
									cplex.addGe(FS[h][g][i][k][p], subexpr3);
								}
							}
						}
					}
				}
			}
			for(int p : this.Ps) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				for(int i : this.N) {
					for(int k : this.K) {
							int [][] A = this.change_v(this.V, i, k);
							for(int h = 0; h < A.length; h ++) {
								for(int g = 0; g < A[h].length; g ++) {
									if(A[h][g] != 0) {
										subexpr1 = this.cplex.sum(subexpr1, this.FS[h][g][i][k][p]);
									}
								}
							}
						}
					}
				cplex.addEq(subexpr1, this.OP1[p]);
			}
			
			//---constraints(1l)---//
			for(int p : this.Pm) {
				for(int k : this.K) {
					for(int h: this.K_0) {
						if(h != k) {
							IloNumExpr subexpr1 = this.cplex.numExpr();
							IloNumExpr subexpr2 = this.cplex.numExpr();
							IloNumExpr subexpr3 = this.cplex.numExpr();
							subexpr1 = this.cplex.prod(this.M, this.XM[h][k][p]);
							subexpr2 = this.cplex.sum(this.tm[h][k], this.Y[k][p]);
							subexpr3 = this.cplex.diff(subexpr2, this.cplex.prod(this.M, this.cplex.diff(1, this.XM[h][k][p])));
							cplex.addLe(FM[h][k][p], subexpr1);
							cplex.addLe(FM[h][k][p], subexpr2);
							cplex.addGe(FM[h][k][p], subexpr3);
						}
					}
				}
			}
			for(int p : this.Pm) {
				IloNumExpr subexpr1 = this.cplex.numExpr();
				for(int k : this.K) {
					for(int h: this.K_0) {
						if(h != k) {
							subexpr1 = this.cplex.sum(subexpr1, this.FM[h][k][p]);
						}
					}
				}
				cplex.addEq(subexpr1, this.OP2[p]);
			}

			//---constraints(1m)---//
			for(int p : this.Ps) {
				cplex.addLe(this.OP1[p], this.Z);
			}
			for(int p : this.Pm) {
				cplex.addLe(this.OP2[p], this.Z);
			}
			
			//---constraints(1.1)---//
			for(int i : this.N) {
				for(int k : this.K) {
					for(int p : this.Ps) {
						IloNumExpr subexpr1 = this.cplex.numExpr();
						IloNumExpr subexpr2 = this.cplex.numExpr();
						IloNumExpr subexpr3 = this.cplex.numExpr();
						IloNumExpr subexpr4 = this.cplex.numExpr();
						int [][] A = this.change_v(this.V, i, k);
						for(int h = 0; h < A.length; h ++) {
							for(int g = 0; g < A[h].length; g ++) {
								if(A[h][g] != 0) {
									subexpr1 = this.cplex.sum(subexpr1, this.XS[h][g][i][k][p]);
								}
							}
						}
						subexpr2 = this.cplex.sum(1, this.cplex.prod(this.M, this.cplex.diff(0, this.RS[i][k][p])));
						subexpr3 = this.cplex.prod(this.M, this.cplex.diff(1, this.RS[i][k][p]));
						subexpr4 = this.cplex.diff(1, this.RS[i][k][p]);
						cplex.addGe(this.Delta1[i][k][p], subexpr2);
						cplex.addLe(this.Delta1[i][k][p], subexpr3);
						cplex.addEq(subexpr1, subexpr4);

					}
				}
			}
			
			//---constraints(1.2)---//
			for(int k : this.K) {
				for(int p : this.Pm) {
					IloNumExpr subexpr1 = this.cplex.numExpr();
					IloNumExpr subexpr2 = this.cplex.numExpr();
					IloNumExpr subexpr3 = this.cplex.numExpr();
					IloNumExpr subexpr4 = this.cplex.numExpr();
					IloNumExpr subexpr5 = this.cplex.numExpr();
					for(int i : this.N) {
						subexpr5 = this.cplex.sum(subexpr5, this.Delta2[i][k][p]);
					}
					for(int h : this.K_0) {
						if(h != k) {
							subexpr1 = this.cplex.sum(subexpr1, this.XM[h][k][p]);
						}
					}
					subexpr2 = this.cplex.sum(1, this.cplex.prod(this.M, this.cplex.diff(0, this.RM[k][p])));
					subexpr3 = this.cplex.prod(this.M, this.cplex.diff(1, this.RM[k][p]));
					subexpr4 = this.cplex.diff(1, this.RM[k][p]);
					cplex.addGe(subexpr5, subexpr2);
					cplex.addLe(subexpr5, subexpr3);
					cplex.addEq(subexpr1, subexpr4);
				}
			}

			cplex.solve();
			cplex.writeSolution("MP_1.txt");
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


				MBT model = new MBT();
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
		        fileWriter = new FileWriter("results\\MBT_"+ conf.sin_mac +"_"+conf.mul_mac +".csv");
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
