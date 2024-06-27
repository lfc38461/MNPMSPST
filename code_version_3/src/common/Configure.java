package common;

/*******
 * This section is define the configure file;
 * The instances and the number of machines can be modified in configure.txt;
********/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.*;
import java.util.*;
public class Configure {
	
	//----------------------------------------------------------
	//--------configure parameters-------------------------------------
	//----------------------------------------------------------
	public ArrayList<String> instance;
	public String sin_mac;
	public String mul_mac;
	
	//initialization;
	public Configure(){
		instance = new ArrayList<String>();
	}
	//read in the configure.txt file;
	public Configure read_configure(String path) throws IOException{
		Scanner cin = new Scanner(new BufferedReader(new FileReader(path)));
		for(int i = 0; i < 4; i++)
			cin.next();
		while(true){
			String str = cin.next();
			if(str.equals("Singlemachine") == true)
				break;
			this.instance.add(str);
		}
		for(int i = 0; i < 1; i++)
			cin.next();
		this.sin_mac = cin.next(); 
		for(int i = 0; i < 2; i++)
			cin.next();
		this.mul_mac = cin.next(); 
		cin.close();
		return this;
	}
}
