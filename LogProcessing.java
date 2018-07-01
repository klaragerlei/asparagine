package Asparagine;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import jbcl.calc.structural.properties.TorsionalAngle;
import jbcl.calc.structural.transformations.ZMatrixToCartesian;
import jbcl.data.formats.PDB;
import jbcl.data.types.Atom;
import jbcl.data.types.PdbAtom;
import jbcl.data.types.Vector3D;

public class LogProcessing {

	private static final double valt = 180 / Math.PI;

	public static void main(String[] args) throws Exception {
		long start = System.currentTimeMillis();

		if (args.length != 1) {
			System.out
					.print("Parameterben meg kellene adni a log file nevit, teljes útvonallal, ahol a log file van, he!");
			System.exit(1);
		}
		String inputFile = args[0]; // The log file

		try {
			File file = new File(inputFile);
			ArrayList<Asparagine> asnlist = processFile(file.getCanonicalPath());
			FileWriter fstream = new FileWriter(file.getAbsolutePath()
					+ ".out.csv");
			BufferedWriter out = new BufferedWriter(fstream);
			printHeadline(out);
			for (Asparagine asn : asnlist) {
				printDihedrals(out, asn);
				/*
				  System.out.println("asn.moleculeName"); for (Atom a :
				  asn.getAtoms()) { System.out.println(a.atomType + " " +
				  a.getId() + " " + a.annotation); }
				 */
			}
			out.close();
			System.out.print(file.getCanonicalPath() + " elkeszult "
					+ (System.currentTimeMillis() - start) + "ms alatt");
		} catch (IOException e1) {
			System.out.print("iras / olvasas hiba. Bocs.\n");
			e1.printStackTrace();
			System.exit(1);
		}
	}

	private static ArrayList<Asparagine> processFile(String filename)
			throws Exception {
		ArrayList<Asparagine> asnlist = new ArrayList<Asparagine>();
		LogReader lr = new LogReader();
		lr.readHeader(filename);
		try {
			ArrayList<double[]> ra = new ArrayList<double[]>();
			ArrayList<double[]> aa = new ArrayList<double[]>();
			ArrayList<double[]> da = new ArrayList<double[]>();
			ArrayList<String> atomnum = new ArrayList<String>();
			ArrayList<String> energy = new ArrayList<String>();
			int molRead = lr.readFromLog(ra, aa, da, atomnum, energy);
			for (int j = 0; j < molRead; j++) {
				Asparagine asn = new Asparagine("Asn" + atomnum.get(j),
						lr.getSd(), ra.get(j), aa.get(j), da.get(j),
						Double.parseDouble(energy.get(j)));
				asnlist.add(asn);
			}
		} catch (FileNotFoundException e) {
			System.out.print(filename + " nincs meg. Pedig kerestem.");
			System.exit(1);
		} catch (IOException e) {
			System.out.print(filename + " beolvasási hiba. Mittudomén.");
			System.exit(1);
		}
		return asnlist;
	}

	private static void printHeadline(BufferedWriter out) throws IOException {
		out.write("Molecule Name;");
		out.write("energy;");
		out.write("omega0;");
		out.write("phi;");
		out.write("psi;");
		out.write("omega;");
		out.write("k1;");
		out.write("k2\n");

	}

	private static void printDihedrals(BufferedWriter out, Asparagine asn)
			throws IOException {
		out.write(asn.moleculeName + ";");
		out.write(asn.getEnergy() + ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.C1, asn.C2, asn.N1, asn.C3)
				+ ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.C2, asn.N1, asn.C3, asn.C4)
				+ ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.N1, asn.C3, asn.C4, asn.N2)
				+ ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.C3, asn.C4, asn.N2, asn.C5)
				+ ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.N1, asn.C3, asn.C6, asn.C7)
				+ ";");
		out.write(valt
				* TorsionalAngle.calculateValue(asn.C3, asn.C6, asn.C7, asn.N3)
				+ "\n");
	}
}
