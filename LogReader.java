package Asparagine;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class LogReader {

	private final StructDescriptor sd = new StructDescriptor();
	private String filename;
	
	private String[] values;

	public StructDescriptor getSd() {
		return sd;
	}

	public void readHeader(String fn) {
		filename = fn;
		FileReader fr = null;
		try {
			fr = new FileReader(filename);
			BufferedReader reader = new BufferedReader(fr);
			String line = null;
			line = reader.readLine(); // keressük az elejit
			while (!line.startsWith(" Symbolic Z-matrix:")) {
				line = reader.readLine();
			}
			line = reader.readLine(); // no még egy sor és kezdödik
			String trimline = reader.readLine().trim();
			while (!trimline.isEmpty()) {
				String an = trimline.split(" ", 2)[0];
				sd.addAtomType(an);
				trimline = reader.readLine().trim();
			}
			line = reader.readLine(); // keressük a struktura leírás elejit
			while (!line.contains("Initial Parameters")) {
				line = reader.readLine();
			}
			line = reader.readLine(); // ugorgyunk
			line = reader.readLine();
			line = reader.readLine();
			line = reader.readLine();
			trimline = reader.readLine().trim();
			while (!trimline.startsWith("---")) {
				String[] tokenized = trimline.split(" ++", 4);
				String[] values = tokenized[2].substring(2).split("[,)]");
				if (tokenized[1].startsWith("R")) {
					sd.addBond(Integer.parseInt(values[0]),
							Integer.parseInt(values[1]));
				}
				if (tokenized[1].startsWith("A")) {
					sd.addAngle(Integer.parseInt(values[0]),
							Integer.parseInt(values[1]),
							Integer.parseInt(values[2]));
				}
				if (tokenized[1].startsWith("D")) {
					sd.addDihAngle(Integer.parseInt(values[0]),
							Integer.parseInt(values[1]),
							Integer.parseInt(values[2]),
							Integer.parseInt(values[3]));
				}
				trimline = reader.readLine().trim();
			}
			fr.close();
		} catch (FileNotFoundException e) {
			System.out.print(filename + " nincs meg. Pedig kerestem.");
			System.exit(1);
		} catch (IOException e) {
			System.out.print(filename + " beolvasási hiba. Mittudomén.");
			System.exit(1);
		}
	}

	public int readFromLog(ArrayList<double[]> distanceList, ArrayList<double[]> angleList,
			ArrayList<double[]> dihAngleList, ArrayList<String> molIdList, ArrayList<String> energyList) throws IOException {
		int nBond = sd.getBonds().size();
		int nAngle = sd.getAngles().size();
		int nDih = sd.getDihAngles().size();
		double[][] ra;
		double[][] aa;
		double[][] da;
		FileReader fr = null;
		fr = new FileReader(filename);
		BufferedReader reader = new BufferedReader(fr);
		String line = null;
		line = reader.readLine(); // keressük a vége elejit
		while (!line.startsWith(" Summary of Optimized")) {
			line = reader.readLine();
		}
		String trimline = reader.readLine().trim();
		int moleculesRead = 0;
		while (!trimline.startsWith("Input orientation:")) {
			String[] atomnum = trimline.split(" ++");
			for (int i = 0; i < atomnum.length; i++) {
				molIdList.add(atomnum[i]);
			}
			moleculesRead += atomnum.length;
			trimline = reader.readLine().trim(); // Eigenvalues
			trimline = trimline.replaceAll("--", " ");
			trimline = trimline.replaceAll("-", " -");
			String[]energy = trimline.split(" ++");
			for (int i = 1; i < energy.length; i++) {
				energyList.add(energy[i]);
			}
			values = reader.readLine().trim().replaceAll("-", " -")
					.split(" ++");
			ra = new double[atomnum.length][nBond];
			aa = new double[atomnum.length][nAngle];
			da = new double[atomnum.length][nDih];
			for (int i = 0; i < nBond; i++) {
				for (int j = 0; j < atomnum.length; j++) {
					ra[j][i] = Double.parseDouble(values[j + 1]);
				}
				values = reader.readLine().trim().replaceAll("-", " -")
						.split(" ++");
			}
			for (int i = 0; i < ra.length; i++) {
				distanceList.add(ra[i]);
			}
			for (int i = 0; i < nAngle; i++) {
				for (int j = 0; j < atomnum.length; j++) {
					aa[j][i] = Double.parseDouble(values[j + 1]);
				}
				values = reader.readLine().trim().replaceAll("-", " -")
						.split(" ++");
			}
			for (int i = 0; i < aa.length; i++) {
				angleList.add(aa[i]);
			}
			for (int i = 0; i < nDih; i++) {
				for (int j = 0; j < atomnum.length; j++) {
					da[j][i] = Double.parseDouble(values[j + 1]);
				}
				trimline = reader.readLine().trim();
				values = trimline.replaceAll("-", " -").split(" ++");
			}
			for (int i = 0; i < ra.length; i++) {
				dihAngleList.add(da[i]);
			}
		}
		return moleculesRead;
	}

}
