package benders;

/*******
 * This section is Run only for LBBD algorithm.
********/


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import ilog.cplex.IloCplex;
import common.Configure;
import common.Read_In;
import java.util.Arrays;
import java.io.File;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import ilog.concert.*;

//the main function;
public class Run {
	
	public static String readPath = "data\\";
	public static void main(String[] args) throws IOException {
		//The instances and the number of machines can be modified in configure.txt;
		Configure conf = new Configure();
		conf.read_configure("configure\\configure.txt");
		//record the computed results by a CSV file;
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
		//Execute each instance in the configuration file;
		for(int i = 0; i < conf.instance.size(); i ++) {
			BMP model = new BMP(conf.instance.get(i));
			model.m1 = Integer.parseInt(conf.sin_mac);
			model.m2 = Integer.parseInt(conf.mul_mac);
			Read_In read_in = new Read_In();
			read_in.read(readPath + conf.instance.get(i) + ".dta", model.ord, model.ARR); 
			System.out.println("Instance : " + conf.instance.get(i) + " M1 = " + model.m1+ " M2 = " + model.m2);
			model.initial(read_in);
			model.model_constract(read_in);
			long starttime = System.currentTimeMillis();
			double bnd = 0.0;
			try {
				LBBD bend = new LBBD(model);
				bend.algorithm();
				bnd = bend.UB;
	
			}catch (IloException e) {
				System.err.println("Concert exception caught '" + e + "' caught");
			}
			long endtime = System.currentTimeMillis();
			System.out.println((endtime-starttime)/1000.00000 + "s");
			model.run_time = (endtime-starttime)/1000.00000;
			model.ARR.clear();
			ArrayList<String>record1 = new ArrayList<String>();
			record1.add(String.valueOf(read_in.K));
			record1.add(String.valueOf(4));
			record1.add(String.valueOf(read_in.num));
			record1.add(String.valueOf(model.m1));
			record1.add(String.valueOf(model.m2));
			record1.add(String.valueOf(model.makespan));
			record1.add(String.valueOf(bnd));
			record1.add(String.valueOf((bnd-model.makespan)/model.makespan));
			record1.add(String.valueOf(model.run_time) + "s");
			record.add(record1);
	

			
		}
		
		//output the computed results to CSV file;
        FileWriter fileWriter = null;
        CSVPrinter csvPrinter = null;
        CSVFormat csvFormat = CSVFormat.DEFAULT;
		try {
		        fileWriter = new FileWriter("results\\LBBD_"+ conf.sin_mac +"_"+conf.mul_mac +".csv");
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
