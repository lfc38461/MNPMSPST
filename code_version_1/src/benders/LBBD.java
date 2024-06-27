package benders;

/*******
 * This section is the code of LBBD.
********/

import java.io.IOException;
import java.util.ArrayList;
import ilog.cplex.IloCplex;
import common.Read_In;
import java.util.Arrays;
import ilog.concert.*;
public class LBBD {
	//parameters;
	public double UB;
	public double LB;
	public BMP mp;
	public ArrayList<BSP> sps;
	public double spub;
	
	//initialization;
	public LBBD(BMP mp) {
		this.UB = Double.MAX_VALUE;
		this.LB = 0.0;
		this.mp = mp;
		this.sps = new ArrayList<BSP>();
	}
	//the LBBD frame;
	public void algorithm() throws IloException{
		try {
			int tempf = 0;
			while((UB - LB > 1e-7)) {
				this.mp.model_solver();
				this.LB = Math.max(LB, this.mp.mpobj);
				BSP sp = new BSP(this.mp);
				sp.solve();
				spub = Math.min(sp.U_B, spub);
				//if the BSP is infeasible, add infeasibility cuts;
				if(sp.obj < 0 || sp.obj > 10000) {
					//single machine;
					for(int p = 0; p < sp.ps; p ++) {
							if(sp.cmaxs.get(p) == sp.obj) {
								IloNumExpr subexpr1 = this.mp.cplex.numExpr();
								for(int i = 0; i < sp.n_s.get(p).size(); i ++) {
									int n = sp.n_s.get(p).get(i);
									int k = sp.k_s.get(p).get(i);
									subexpr1 = this.mp.cplex.sum(subexpr1, 
											   this.mp.cplex.abs(
													   this.mp.cplex.diff(this.mp.Delta1[n][k][p], sp.delt1[n][k][p])));
									this.mp.cplex.addLe(this.mp.Delta1[n][k][p], sp.delt1[n][k][p]);//
								}
								this.mp.cplex.addGe(subexpr1, 1 );
							}	
					}

					//multi machine;
					for(int p = 0; p < sp.pm; p ++) {
						int cnt = 0;
						ArrayList<Integer> n_store = new ArrayList<Integer>();
						ArrayList<Integer> k_store = new ArrayList<Integer>();
						for(int k = 0; k < sp.K; k ++) {
							if(sp.y[k][p] > 0) {
								for(int n = 0; n < sp.N; n ++) {
									if(sp.delt2[n][k][p] * this.mp.wm[n][k] == sp.y[k][p]) {
										cnt += sp.delt2[n][k][p];
										n_store.add(n);
										k_store.add(k);
										break;
									}
								}
							}
						}
						if(sp._cmaxm[p] == sp.obj) { 
							IloNumExpr subexpr1 = this.mp.cplex.numExpr();
							for(int i = 0; i < n_store.size(); i ++) {
								int a = n_store.get(i);
								int b = k_store.get(i);
								subexpr1 = this.mp.cplex.sum(subexpr1, 
										   this.mp.cplex.abs( 
												   this.mp.cplex.diff(this.mp.Delta2[a][b][p], sp.delt2[a][b][p])));
								this.mp.cplex.addLe(this.mp.Delta2[a][b][p], sp.delt2[a][b][p]);//
							}
							this.mp.cplex.addGe(subexpr1, 1);
						}
					}	
				}
				else {
					//if the BSP is feasible, add optimality cuts;
					this.UB = Math.min(this.UB, sp.obj);
					//single machine;
					for(int p = 0; p < sp.n_s.size(); p ++) {
						IloNumExpr subexpr1 = this.mp.cplex.numExpr();
						IloNumExpr subexpr2 = this.mp.cplex.numExpr();
						IloNumExpr subexpr3 = this.mp.cplex.numExpr();
						for(int n = 0; n < sp.n_s.get(p).size(); n ++) {
							int nn = sp.n_s.get(p).get(n);
							int kk = sp.k_s.get(p).get(n);
							subexpr1 = this.mp.cplex.sum(subexpr1, 
											this.mp.cplex.prod(	this.mp.ws[nn][kk],
													this.mp.cplex.diff(sp.delt1[nn][kk][p],
													 this.mp.Delta1[nn][kk][p])));	
						}
						for(int n = 0; n < sp.n_s.get(p).size(); n ++) {
							int nn = sp.n_s.get(p).get(n);
							int kk = sp.k_s.get(p).get(n);
							subexpr2 = this.mp.cplex.sum(subexpr2,
										   	this.mp.cplex.prod(sp.ES.get(p).get(n),
										   		this.mp.cplex.diff(sp.xs[nn][kk][p], this.mp.XS[nn][kk][p])));
						}
						subexpr3 = this.mp.cplex.diff(sp._cmaxs[p], 
										this.mp.cplex.sum(subexpr1, subexpr2));
						this.mp.cplex.addGe(this.mp.C, subexpr3);
					}
					
					
					//multi machine;
					for(int p = 0; p < sp.k_m.size(); p ++) {
						IloNumExpr subexpr1 = this.mp.cplex.numExpr();
						IloNumExpr subexpr2 = this.mp.cplex.numExpr();
						IloNumExpr subexpr3 = this.mp.cplex.numExpr();

						for(int k = 0; k < sp.k_m.get(p).size(); k ++) {
							int kk = sp.k_m.get(p).get(k);
							subexpr1 = this.mp.cplex.sum(subexpr1,
									   		this.mp.cplex.diff(sp.y[kk][p], this.mp.Y[kk][p]));
						}
						for(int k = 0; k < sp.K; k ++) {
							if(sp.xm[k][p] > 0.0) {
								subexpr2 = this.mp.cplex.sum(subexpr2,
										this.mp.cplex.prod(sp.em_b[k],
									   		this.mp.cplex.diff(sp.xm[k][p], this.mp.XM[k][p])));
							}else {
								subexpr2 = this.mp.cplex.sum(subexpr2,
										this.mp.cplex.prod(sp.em_s[k],
									   		this.mp.cplex.diff(sp.xm[k][p], this.mp.XM[k][p])));
							}
						}
						subexpr3 = this.mp.cplex.diff(sp._cmaxm[p], 		
										this.mp.cplex.sum(subexpr1, subexpr2));
						this.mp.cplex.addGe(this.mp.C, subexpr3);						
					}
				}

				tempf ++;
				if(tempf > 0) {
				}
				this.mp.cplex = sp.mp.cplex;
			}
			System.out.println("LB = " + LB);
			System.out.println("UB = " + UB);
		}
			catch(IOException e) {
				System.err.println("Concert exception caught '" + e + "' caught");
			}
			
		
	}
}
