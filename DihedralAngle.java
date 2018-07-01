
package Asparagine;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import jbcl.calc.structural.properties.TorsionalAngle;

	public class DihedralAngle {

		static BufferedWriter out;
		private static final double valt = 180/Math.PI;

		public static void main(String[] args) {
			long start = System.currentTimeMillis();

			if (args.length!=1){
				System.out.print("Parameterben meg kellene adni a mappa nevit, ahol a pdb-k vannak, he!");
				System.exit(1);
			}
			String inputPDBFolder = args[0]; // A folder to read pdb files

			File folder = new File(inputPDBFolder);
			File[] listOfFiles = folder.listFiles();

			try {
				FileWriter fstream = new FileWriter(inputPDBFolder + "\\out.txt");
				out = new BufferedWriter(fstream);

				for (File file : listOfFiles)
					if (!file.isDirectory() && file.getName().endsWith(".pdb"))
						processFile(file.getCanonicalPath());
				out.close();
			} catch (IOException e1) {
				System.out.print("iras / olvasas hiba. Bocs.\n");
				e1.printStackTrace();
				System.exit(1);
			}
			System.out.print(inputPDBFolder + "\\out.txt elkeszult "+ (System.currentTimeMillis()-start) + "ms alatt");
		}

		private static void processFile(String filename) throws IOException {
			Asparagine asn = new Asparagine(filename);
			out.write(filename + "\n");
			out.write(asn.moleculeName + "\n");
			out.write("dihedral omega0="
					+ valt*TorsionalAngle.calculateValue(asn.C1, asn.C2, asn.N1, asn.C3)
					+ "\n");
			out.write("dihedral phi   ="
					+ valt*TorsionalAngle.calculateValue(asn.C2, asn.N1, asn.C3, asn.C4)
					+ "\n");
			out.write("dihedral psi   ="
					+ valt*TorsionalAngle.calculateValue(asn.N1, asn.C3, asn.C4, asn.N2)
					+ "\n");
			out.write("dihedral omega ="
					+ valt*TorsionalAngle.calculateValue(asn.C3, asn.C4, asn.N2, asn.C5)
					+ "\n");
	/*		out.write("dihedral C2,N1,C3,C6="
					+ valt*TorsionalAngle.calculateValue(asn.C2, asn.N1, asn.C3, asn.C6)
					+ "\n");
	*/
			out.write("dihedral k1    ="
					+ valt*TorsionalAngle.calculateValue(asn.N1, asn.C3, asn.C6, asn.C7)
					+ "\n");
			out.write("dihedral k2    ="
					+ valt*TorsionalAngle.calculateValue(asn.C3, asn.C6, asn.C7, asn.N3)
					+ "\n");

		}

	}
