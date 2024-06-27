package benders;

/*******
 * This section is the code of Benders SubProblem.
********/

import java.util.ArrayList;
import java.util.Collections;
import ilog.cplex.IloCplex;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;

public class BSP {
	//parameters
	public double obj; 
	public double U_B;
	public Boolean infeasible;
	public BMP mp;
	public int N;
	public int K;
	public int ps;
	public int pm;
	public double xs [][][];
	double xm [][];
	double delt1[][][];
	double delt2[][][];
	double y    [][];
	double _cmaxs[];
	double _cmaxm[];
	double yk [][][];
	int em_b[];
	int em_s[];
	ArrayList<Integer> cont_delt;      //record the sum of delta on each machine;
	ArrayList<Double> cmaxs;
	ArrayList<Double> cmaxm;
	ArrayList<ArrayList<Integer>> n_s;  //the i coordinate of each single machine with x=1;
	ArrayList<ArrayList<Integer>> k_s;	//the k coordinate of each single machine with x=1;
	ArrayList<ArrayList<Integer>> ES;	//max setup time between every 2 type on single machine;
	ArrayList<ArrayList<Integer>> ES_s; //min setup time between every 2 type on single machine;
	ArrayList<ArrayList<Integer>> k_m;
	ArrayList<ArrayList<Integer>> EM;
	ArrayList<ArrayList<Integer>> EM_s;
	
	
	
	
	//initialization;
	public BSP(BMP mp) {
		try {
			this.mp = mp;
			this.N = this.mp.N_0.size();
			this.K = this.mp.K_0.size();
			this.ps = this.mp.Ps.size();
			this.pm = this.mp.Pm.size();
			this.xs = new double[N][K][ps];
			this.xm = new double[K][pm];
			this.delt1 = new double [N][K][ps];
			this.delt2 = new double [N][K][pm];
			this.yk = new double [N][K][pm];
			this.y     = new double [K][pm];
			this._cmaxs = new double[ps];
			this._cmaxm = new double[pm];
			this.cmaxs = new ArrayList<Double>();
			this.cmaxm = new ArrayList<Double>();
			this.cont_delt = new ArrayList<Integer>();
			this.U_B = Double.MAX_VALUE;
			n_s = new ArrayList<ArrayList<Integer>>();
			k_s = new ArrayList<ArrayList<Integer>>();
			EM  = new ArrayList<ArrayList<Integer>>();
			EM_s  = new ArrayList<ArrayList<Integer>>();
			k_m = new ArrayList<ArrayList<Integer>>();
			ES  = new ArrayList<ArrayList<Integer>>();
			ES_s  = new ArrayList<ArrayList<Integer>>();
			
			for(int i = 0; i < N; i ++) {
				for(int j = 0; j < K; j ++) {
					for(int k = 0; k < ps; k++) {
						xs[i][j][k] = this.mp.cplex.getValue(this.mp.XS[i][j][k]);
						delt1[i][j][k] = this.mp.cplex.getValue(this.mp.Delta1[i][j][k]);
					}
				}
			}
			for(int i = 0; i < N; i ++) {
				for(int j = 0; j < K; j ++) {
					for(int k = 0; k < pm; k++) {
						delt2[i][j][k] = this.mp.cplex.getValue(this.mp.Delta2[i][j][k]);
						yk[i][j][k]	   = delt2[i][j][k] * mp.wm[i][j];
					}
				}
			}
			for(int i = 0; i < K; i ++) {
				for(int j = 0; j < pm; j ++) {
					xm[i][j] = this.mp.cplex.getValue(this.mp.XM[i][j]);
					y[i][j]  = this.mp.cplex.getValue(this.mp.Y[i][j]);
				}
			}
			
		
		}catch(IloException e) {
			System.err.println("Concert exception caught '" + e + "' caught");
		}
	}
	public void set_mp(BMP mp) {
		this.mp = mp;
	}
	
	
	//function for solve ATSP in subproblem;
	public double atsp(ArrayList<Integer> org_set0, ArrayList<Integer> set_k, ArrayList<Integer> set_k0) {
		IloCplex SP_cplex;
		IloObjective SP_obj;
		IloNumVar [][]   SP_X;
		IloNumVar []     SP_Y;
		double sp_cost = 0;
		try {
			SP_cplex = new IloCplex();
			SP_cplex.setParam(IloCplex.Param.Simplex.Tolerances.Optimality, 1e-9);
			SP_cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-9);
			SP_cplex.setParam(IloCplex.DoubleParam.TimeLimit, 50);
			SP_cplex.setOut(null);
			
			SP_X = new IloNumVar[set_k0.size()][set_k0.size()];
			for(int i = 0; i < set_k0.size(); i ++) {
				for(int j = 0; j < set_k0.size(); j ++) {
					SP_X[i][j] = SP_cplex.numVar(0.0, 1, IloNumVarType.Int, "SP_X["+i+"]["+j+"]");
				}
			}
			
			SP_Y = new IloNumVar[set_k0.size()];
			for(int i = 0; i < set_k0.size(); i ++) {
				SP_Y[i]  = SP_cplex.numVar(0.0, set_k.size(), IloNumVarType.Int, "SP_Y["+i+"]");
			}
			
			IloNumExpr objexpr = SP_cplex.numExpr();
			for(int i : set_k0) {
				for(int k :set_k0) {
					objexpr = SP_cplex.sum(objexpr, SP_cplex.prod(SP_X[i][k], this.mp.tm[org_set0.get(i)][org_set0.get(k)]));
				}
			}
			SP_obj = SP_cplex.addMinimize(objexpr);

			
			//cons1
			IloNumExpr subexpr11 = SP_cplex.numExpr();
			IloNumExpr subexpr22 = SP_cplex.numExpr();
			for(int i : set_k) {
				subexpr11 = SP_cplex.sum(subexpr11, SP_X[0][i]);
				subexpr22 = SP_cplex.sum(subexpr22, SP_X[i][0]);
				SP_cplex.addEq(SP_X[i][i], 0);
			}
			SP_cplex.addEq(subexpr11, 1); SP_cplex.addEq(subexpr22, 1);
			
			//cons2
			for(int i : set_k) {
				IloNumExpr subexpr1 = SP_cplex.numExpr();
				IloNumExpr subexpr2 = SP_cplex.numExpr();
				for(int j : set_k0) {
					subexpr1 = SP_cplex.sum(subexpr1, SP_X[i][j]);
				}
				for(int j : set_k0) {
					subexpr2 = SP_cplex.sum(subexpr2, SP_X[j][i]);
				}
				SP_cplex.addEq(subexpr1, 1); SP_cplex.addEq(subexpr2, 1);
			}
			
			//cons3
			for(int i : set_k) {
				for(int j : set_k) {
					if(i != j) {
						IloNumExpr subexpr1 = SP_cplex.numExpr();
						subexpr1 = SP_cplex.sum(
								SP_cplex.diff(SP_Y[i], SP_Y[j]), SP_cplex.prod(set_k0.size(), SP_X[i][j])
								);
						SP_cplex.addLe(subexpr1, set_k0.size() - 1);
					}
				}
			}
			
			SP_cplex.solve();
			sp_cost = SP_cplex.getObjValue();
			SP_cplex.writeSolution("SP_MULTI.txt");
			SP_cplex.end();
			}catch(IloException e) {
			System.err.println("Concert exception caught '" + e + "' caught");
			}
		return sp_cost;
	}
	
	
	
	
	//function for solve subproblem;
	public void solve() {
		//solve for single machine;
		for(int p = 0; p < this.ps; p ++) {
			ArrayList<Integer> n_cr = new ArrayList<Integer>();
			ArrayList<Integer> k_cr = new ArrayList<Integer>();
			ArrayList<Double> cost = new ArrayList<Double>();			
			int cnt = 0;
			for(int n = 0; n < this.N; n ++) {
				for(int k = 0; k < this.K; k ++) {
					if(Math.abs(this.xs[n][k][p]-1.0) <= 0.00001) {
						n_cr.add(n);
						k_cr.add(k);
						cnt += this.delt1[n][k][p];
					}
				}
			}
			ArrayList<Integer> biggest_t = new ArrayList<Integer>();
			ArrayList<Integer> smallest_t = new ArrayList<Integer>();
			for(int i = 0; i < n_cr.size(); i ++) {
				int a = n_cr.get(i);
				int b = k_cr.get(i);
				int temp = 0; 
				int temp1 = Integer.MAX_VALUE;
				for(int ii = 1; ii < this.N; ii ++) {
					for(int jj = 1; jj < this.K; jj ++) {
						if(this.mp.ts[ii][jj][a][b] > temp) {
							temp = this.mp.ts[ii][jj][a][b];
						}
						if(this.mp.ts[ii][jj][a][b] < temp1) {
							temp1 = this.mp.ts[ii][jj][a][b];
						}
					}
				}
				if(temp > this.mp.ts[0][0][a][b]) {
					biggest_t.add(temp);
				}
				else {
					biggest_t.add(this.mp.ts[0][0][a][b]);
				}
				if(temp1 < this.mp.ts[0][0][a][b]) {
					smallest_t.add(temp1);
				}
				else {
					smallest_t.add(this.mp.ts[0][0][a][b]);
				}
			}
			for(int i = 0; i < n_cr.size(); i ++) {
				Double count = 0.0;
				for(int j = 0; j < n_cr.size(); j ++) {
					if(i == j) {
						count +=  mp.ts[0][0][n_cr.get(j)][k_cr.get(j)] 
								+ mp.ws[n_cr.get(j)][k_cr.get(j)] * this.delt1[n_cr.get(j)][k_cr.get(j)][p];
					}
					else {
						int temp = 0;
						for(int ii = 1; ii < this.N ; ii ++) {
							for(int jj = 1; jj < this.K ; jj ++)
								if(mp.ts[ii][jj][n_cr.get(j)][k_cr.get(j)] > 0) {
									temp = mp.ts[ii][jj][n_cr.get(j)][k_cr.get(j)] ;
									break;
								}
						}
						count += temp
								+ mp.ws[n_cr.get(j)][k_cr.get(j)] * this.delt1[n_cr.get(j)][k_cr.get(j)][p];
					}
				}
				cost.add(count);
			}

			Collections.sort(cost);
			
			if(n_cr.size() != 0) {
				this.cont_delt.add(cnt);
				this._cmaxs[p] = cost.get(0);
				this.cmaxs.add(cost.get(0));
				this.n_s.add(n_cr);
				this.k_s.add(k_cr);
				this.ES.add(biggest_t);
				this.ES_s.add(smallest_t);
			}
			else {
				this._cmaxs[p] = -1;
			}

		}
		
		//solve for multi machine;
		for(int p = 0; p < this.pm; p ++) {
			ArrayList<Integer> k_cr = new ArrayList<Integer>();
			ArrayList<Double> cost = new ArrayList<Double>();
			for(int k = 0; k < this.K; k ++) {
				if(Math.abs(this.xm[k][p] - 1) < 0.0001) {
					k_cr.add(k);
				}
			}
			ArrayList<Integer> biggest_t = new ArrayList<Integer>();
			ArrayList<Integer> smallest_t = new ArrayList<Integer>();
			for(int i = 0; i < k_cr.size(); i ++) {
				int a = k_cr.get(i);
				int temp = 0; 
				int temp1 = Integer.MAX_VALUE;
				for(int ii = 1; ii < this.K; ii ++) {
					if(this.mp.tm[ii][a] > temp) {
						temp = this.mp.tm[ii][a];
					}
					if(this.mp.tm[ii][a] < temp1) {
						temp1 = this.mp.tm[ii][a];
					}
				}
				if(temp > this.mp.tm[0][a]) {
					biggest_t.add(temp);
				}
				else {
					biggest_t.add(this.mp.tm[0][a]);
				}
				if(temp1 < this.mp.tm[0][a]) {
					smallest_t.add(temp1);
				}
				else {
					smallest_t.add(this.mp.tm[0][a]);
				}
			}
			double mlticost = 0;
			ArrayList<Integer> set_k = new ArrayList<Integer>();
			ArrayList<Integer> set_k0 = new ArrayList<Integer>();
			ArrayList<Integer> org_set0 = new ArrayList<Integer>();
			set_k0.add(0);  org_set0.add(0);
			for(int i = 0; i < k_cr.size(); i ++) {
				set_k.add(i + 1);
				set_k0.add(i + 1);
				org_set0.add(k_cr.get(i));					
			}
			double time_k = 0;
			for(int j = 0; j < k_cr.size(); j ++) { 
				double temp_k = 0;
				for(int n = 0; n < this.N; n ++){ 
					if(this.delt2[n][k_cr.get(j)][p] > 0) {
						if(this.delt2[n][k_cr.get(j)][p] * this.mp.wm[n][k_cr.get(j)] > temp_k){
							temp_k = this.delt2[n][k_cr.get(j)][p] * this.mp.wm[n][k_cr.get(j)];
						}
					}
					}
				time_k += temp_k;
			}
			
			mlticost = atsp(org_set0, set_k, set_k0);
			double aaacost = mlticost + time_k;
			cost.add(aaacost);

			if(k_cr.size() != 0) {
				this._cmaxm[p] = cost.get(0);
				this.cmaxm.add(cost.get(0));
				this.k_m.add(k_cr);
				this.EM.add(biggest_t);
				this.EM_s.add(smallest_t);
			}
			else {
				this._cmaxm[p] = -1;
			}

		}
		this.obj = 0;
		for(int i = 0; i < this.cmaxs.size(); i ++) {
			if(this.cmaxs.get(i) > this.obj) {
				this.obj = this.cmaxs.get(i);
			}
		}
		for(int i = 0; i < this.cmaxm.size(); i ++) {
			if(this.cmaxm.get(i) > this.obj) {
				this.obj = this.cmaxm.get(i);
			}			
		}
		this.U_B = this.obj;
		this.em_b = new int[mp.tm.length];
		this.em_s = new int[mp.tm.length];
		int row_b_s [][] = new int [mp.tm.length][2];
		int col_b_s [][] = new int [mp.tm.length][2];
		for(int i = 1; i < mp.tm.length; i ++) {
			int temp = 0; 
			int temp1 = Integer.MAX_VALUE;
			int tempc = 0; 
			int tempc1 = Integer.MAX_VALUE;
			for(int j = 1; j < mp.tm[1].length; j ++) {
				if(mp.tm[i][j] > temp) temp = mp.tm[i][j];
				if(mp.tm[i][j] < temp1) temp1 = mp.tm[i][j];
				if(mp.tm[j][i] > temp) tempc = mp.tm[j][i];
				if(mp.tm[j][i] < temp1) tempc1 = mp.tm[j][i];
			}
			row_b_s[i][0] = temp1; row_b_s[i][1] = temp; 
			col_b_s[i][0] = tempc1; col_b_s[i][1] = tempc;
		}
		for(int i = 0; i < em_b.length; i ++) {
			em_b[i] = Math.max(row_b_s[i][1], col_b_s[i][1]);
			em_s[i] = Math.min(row_b_s[i][0], col_b_s[i][0]);
		}	
	}
	
	

	
}
