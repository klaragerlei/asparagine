/**
 * 
 */
package scaninput_generator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

/**
 * @author Klári
 *
 */
public class InputGen {

	/**
	 * @param args
	 */
	private static String filename;
	private static int resolution;
	private static ArrayList<String> dihedralAngles=new ArrayList<String>();
	private static ArrayList<String> lines=new ArrayList<String>();
	private static int dimensions;
	private static int steps;
	private static int numberOfFiles;
	
	public static void main(String[] args) {
		getParameters(args);
		read();
		write();

	}

	private static final void getParameters(String[] args)
	{
		if(args.length<3)
		{
			System.out.println("Usage:java InputGen filename resolution dihedral angles[]"
					+ "(resolution is the degrees the scan takes in one step,"
					+ " dihedral angles should be listed like D12 D6 [D9,D10,D11] etc)");
			System.exit(0);
		}
		filename = args[0];
		resolution = Integer.parseInt(args[1]);
		dimensions = args.length-2;
		steps = 360/resolution;
		for(int i=0;i<(dimensions);i++)
		{
			dihedralAngles.add(i, args[i+2]);
		}
		numberOfFiles = (int) Math.pow(steps,args.length-2);
	}
	
	private static final void read() 
	{
		 BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(filename));

		    try {
		        String line = br.readLine();
		        lines.add(line);

		        while (line != null) {
		            line = br.readLine();
			        lines.add(line);
		            
		        }

		    } finally {
		        br.close();
		    }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	}


	private static final void write()
	{
		PrintWriter writer;
		int[]angle=new int[dimensions];
		try {
            for (int f=0; f<dimensions;f++)
            {
            	angle[f]=0;
            }
	    for(int i = 0; i<numberOfFiles;i++)
	    {
	    	String fileNum="000000000"+new Integer(i).toString();
	    	fileNum=fileNum.substring(fileNum.length()-5);
	    	writer = new PrintWriter(filename+".out\\file"+fileNum+".gjf", "UTF-8");
	    	for(int j = 0; j<lines.size();j++)
	    	{
	    		String line=lines.get(j);
	    		line=processLine(line,angle);
	    		writer.println(line);
	    	}
		writer.close();
        increaseAngles(angle);    
	    }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static String processLine(String line,int[] angles)
	{
		if (line!=null && line.trim()!=null && line.trim().startsWith("D")){
			String[]sp=line.trim().split(" +");
			int i=0;
			boolean found=false;
			while(i<dimensions && !found)
			{		
				if (isToRotate(sp[0],dihedralAngles.get(i)))
				{
					Double a=Double.parseDouble(sp[1])+angles[i];
					while (a>360)
					{
						a-=360;
					}
					line=line.replace(sp[1], a.toString());
					found=true;
				}
				else
				{
					i++;
				}
			}	
		}
		return line;
	}
	
	private static boolean isToRotate(String split, String rotateGroup)
	{
		boolean toRotate = false;
		if(rotateGroup.startsWith("[")) 
		{
			String[]sp=rotateGroup.substring(1, rotateGroup.length()-1).split(",");
			for (int i=0;i<sp.length;i++)
			{
				if(sp[i].equals(split))
				{
					toRotate = true;
		//			break;
				}
			}
		}
		else 
		{
			toRotate = split.equals(rotateGroup);
		}
	    return toRotate;
	}
	
	
	private static void increaseAngles(int[] angle)
	{
		int size=angle.length;
		int i=0;
		boolean done=false;
		while (!done & i<size)
		{
			angle[i]=angle[i]+resolution;
			if (angle[i]<360)
			{
				done=true;
			}
			else
			{
				angle[i]=0;
				i=i+1;
			}
		}
	}
	
}
