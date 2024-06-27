package common;

/*******
 * This section is define the instance-read-in function;
********/

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Read_In {

	//preliminary variable;
	public int    machine1;	//number of single machines;
	public int    machine2;	//number of multi machines;
	public int    N;   		//number of patterns£»
	public int    K;   		//number of type;
	public int    num; 		//total number of jobs£»
	public ArrayList<Integer> typ_num;
	public double duration;
	
	//initialization;
	public Read_In() {
		this.machine1 = 0;
		this.machine2 = 0;
		this.N        = 0;
		this.K		  = 0;
		this.duration = 400;	//unused;
		this.num      = 0;
	}
	 
//read in the instances;	
public void read(String path, int ord, ArrayList<ArrayList<Integer>> arr) throws IOException{
	  File file = new File(path);
	  if (!file.exists()) {
	   System.out.println("File not exist!");
	   System.exit(0);
	  }
	  BufferedReader br = new BufferedReader(new FileReader(path));	  
	  String temp = "";
	  typ_num = new ArrayList<Integer>();
	  while ((temp = br.readLine()) != null) {
		if (temp.startsWith("index,"))
		{
			continue;
		}
		else if(!temp.equals("") )
		{
			String [] str1 = temp.split(" ");
			ArrayList<Integer> array = new ArrayList<Integer>();
			array.add(Integer.parseInt(str1[0]));	 //index;
			array.add(Integer.parseInt(str1[1]));    //processing time;
			array.add(Integer.parseInt(str1[2]));    //class;
			array.add(Integer.parseInt(str1[3]));	 //direction;
			array.add(Integer.parseInt(str1[4]));	 //setup;
			array.add(Integer.parseInt(str1[5]));	 //inisetup;
			array.add(Integer.parseInt(str1[6]));	 //quantity;

			arr.add(array);
			this.num += Integer.parseInt(str1[6]);
			if(this.N < Integer.parseInt(str1[3]))
				N = Integer.parseInt(str1[3]);
			if(this.K <Integer.parseInt(str1[2]))
				K = Integer.parseInt(str1[2]);
  	}
	  }
	  br.close(); 
	  N ++;
	  K ++;

	 }
	}

