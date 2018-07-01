package Asparagine;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.ListIterator;

import Asparagine.StructDescriptor.BondDesc;

import jbcl.chemistry.Molecule;
import jbcl.data.dict.AtomType;
import jbcl.data.dict.BondType;
import jbcl.data.types.Atom;
import jbcl.data.types.Bond;
import jbcl.data.types.Vector3D;
import jbcl.calc.structural.transformations.SphericalToCartesian;
import jbcl.calc.structural.transformations.ZMatrixToCartesian;

public class Asparagine extends Molecule {

	private static final double valt = 180 / Math.PI;

//	private static final double deg2rad = Math.PI / 180;


	private String title = null;
	public Atom C1, C2, C3, C4, C5, C6, C7, N1, N2, N3, O1, O2, O3;
	private double energy;

	// Olvassuk ki egy pdb file-ból
	public Asparagine(String filename) {
		readPdb(filename);
		annotateAtoms();
	}

	// Csináljuk a struktúra leíróból.
	public Asparagine(String name, StructDescriptor sd, double[] r, double[] a, double[] d,
			double energ) throws Exception {
		if (sd == null)
			throw new Exception(
					"Elébb a struktúra leírást kéne feltölteni, mielõtt Asparagint kezdenél szintetizálni");
		int i = 1;
		for (String an : sd.getAtomTypes()) {
			addAtom(new Atom(an, i++, 0, i, i * i));
		}
		for (StructDescriptor.BondDesc bd : sd.getBonds()) {
			optBindAtoms(bd.getA1(), bd.getA2());
		}
		moleculeName = name;
		energy = energ;
		annotateAtoms();
		amendByZmatix(sd, r, a, d);
	}

	public Atom getAtom(int id) {
		ListIterator<Atom> it = getAtoms().listIterator();
		while (it.hasNext()) {
			Atom a = it.next();
			if (a.getId() == id)
				return a;
		}
		return null;
	}

	public double getEnergy() {
		return energy;
	}

	public void setEnergy(double energy) {
		this.energy = energy;
	}

	private void optBindAtoms(int id1, int id2) {
		if (id1 < id2) { // don't try to store bonds both ways
			Atom a1 = getAtom(id1);
			Atom a2 = getAtom(id2);
			if (a1 != null && a2 != null) {
				bindAtoms(a1, a2, BondType.UNKNOWN);
			}
		}
	}

	private Atom getOtherEndOfBond(Atom a, int bondId) {
		if (bondId >= getBonds(a).size())
			return null;
		Bond b = ((Bond) getBonds(a).get(bondId));
		return getOtherEndOfBond(a, b);
	}

	private Atom getOtherEndOfBond(Atom a, Bond b) {
		Atom x = b.secondPartner;
		if (a == x) {
			x = b.firstPartner;
		}
		return x;
	}

	private void readPdb(String filename) {
		FileReader fr = null;
		try {
			fr = new FileReader(filename);
			BufferedReader reader = new BufferedReader(fr);
			String line = null;
			line = reader.readLine(); // TITLE
			if (!line.startsWith("TITLE")) {
				System.out.print("Furcsa, nem TITLE sorral kezdõdik");
				System.exit(1);
			}
			title = line.substring(6).trim();
			line = reader.readLine(); // REMARK
			if (!line.startsWith("REMARK")) {
				System.out.print("Furcsa, nem REMARK a második sor");
				System.exit(1);
			}
			moleculeName = title;
			while ((line = reader.readLine()) != null
					&& !line.startsWith("END")) {
				if (!line.startsWith("HETATM")) {
					System.out
							.print("Furcsa, nem HETATM sor jött, pedig azt vártam");
					System.exit(1);
				}
				// if (line.charAt(77)=='H') continue; // Skip Hydrogenes
				addAtom(new Atom(
						AtomType.valueOf(line.substring(12, 16).trim()),
						new Integer(line.substring(6, 11).trim()).intValue(),
						"",
						new Double(line.substring(30, 38).trim()).doubleValue(),
						new Double(line.substring(38, 46).trim()).doubleValue(),
						new Double(line.substring(46, 54).trim()).doubleValue()));
				// System.out.print(line+"\n");
			}
			while ((line = reader.readLine()) != null) {
				if (!line.startsWith("CONECT")) {
					System.out
							.print("Furcsa, nem CONECT sor jött, pedig azt vártam");
					System.exit(1);
				}
				line = line
						.concat("                                                                  ");
				int myId = new Integer(line.substring(6, 11).trim()).intValue();
				if (getAtom(myId) != null) {
					int neighbour1, neighbour2, neighbour3, neighbour4;
					neighbour1 = new Integer(line.substring(11, 16).trim())
							.intValue();
					String s = line.substring(16, 21).trim();
					if (s.length() != 0) {
						optBindAtoms(myId, neighbour1); // If only one neighbour
														// then its a H and we
														// skip those
						neighbour2 = new Integer(s).intValue();
						optBindAtoms(myId, neighbour2);
						s = line.substring(21, 26).trim();
						if (s.length() != 0) {
							neighbour3 = new Integer(s).intValue();
							optBindAtoms(myId, neighbour3);
							s = line.substring(26, 31).trim();
							if (s.length() != 0) {
								neighbour4 = new Integer(s).intValue();
								optBindAtoms(myId, neighbour4);
							}
						}
					}
				}
				// System.out.print(line+"\n");
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

	private int bondsNoH(Atom a) {
		LinkedList<Bond> lb = getBonds(a);
		int noH = 0;
		for (Bond b : lb) {
			if (getOtherEndOfBond(a, b).atomType != AtomType.H)
				noH++;
		}
		return noH;
	}

	private Atom neighbourOfType(Atom a, AtomType at) {
		int i = 0;
		Atom b = getOtherEndOfBond(a, i++);
		while (b != null && b.atomType != at) {
			b = getOtherEndOfBond(a, i++);
		}
		return b;
	}

	private Atom neighbourOfType(Atom a, AtomType at, Atom excluded) {
		int i = 0;
		Atom b = getOtherEndOfBond(a, i++);
		while (b != null && (b.atomType != at || b == excluded)) {
			b = getOtherEndOfBond(a, i++);
		}
		return b;
	}

	private void annotateAtoms() {
		LinkedList<Atom> list = getAtoms();
		ListIterator<Atom> lit = list.listIterator();
		Atom a = null, b = null, c = null;
		while (lit.hasNext()) {
			a = lit.next();
			if (a.atomType == AtomType.C && bondsNoH(a) == 1
					&& neighbourOfType(a, AtomType.C) != null) {
				a.annotation = "C_1";
				C1 = a;
				break;
			}
		}
		b = neighbourOfType(a, AtomType.C);
		a = b;
		a.annotation = "C_2";
		C2 = a;
		b = neighbourOfType(a, AtomType.O);
		b.annotation = "O_1";
		O1 = b;
		b = neighbourOfType(a, AtomType.N);
		b.annotation = "N_1";
		N1 = b;
		a = b;
		b = neighbourOfType(a, AtomType.C, C2);
		b.annotation = "C_3";
		C3 = b;
		a = b;
		b = neighbourOfType(a, AtomType.C);
		c = neighbourOfType(a, AtomType.C, b);
		if (bondsNoH(b) == 2) {
			b.annotation = "C_6";
			C6 = b;
			c.annotation = "C_4";
			C4 = c;
		} else {
			b.annotation = "C_4";
			C4 = b;
			c.annotation = "C_6";
			C6 = c;
		}
		a = C4;
		b = neighbourOfType(a, AtomType.O);
		b.annotation = "O_2";
		O2 = b;
		b = neighbourOfType(a, AtomType.N);
		b.annotation = "N_2";
		N2 = b;
		a = b;
		b = neighbourOfType(a, AtomType.C, C4);
		b.annotation = "C_5";
		C5 = b;
		a = C6;
		b = neighbourOfType(a, AtomType.C, C3);
		b.annotation = "C_7";
		C7 = b;
		a = C7;
		b = neighbourOfType(a, AtomType.O);
		b.annotation = "O_3";
		O3 = b;
		b = neighbourOfType(a, AtomType.N);
		b.annotation = "N_3";
		N3 = b;
	}

	// =================== LOG file processing =================/


	private void amendByZmatix(StructDescriptor sd, double[] r, double[] a, double[] d) {
		ArrayList<BondDesc> bonds = sd.getBonds();
		if (bonds.size() != r.length) {
			System.out.print("Kötés lista hossza rossz: desc=" + bonds.size()
					+ " r=" + r.length);
			System.exit(1);
		}
		if (sd.getAngles().size() != a.length) {
			System.out.print("Szög lista hossza rossz: desc="
					+ sd.getAngles().size() + " a=" + a.length);
			System.exit(1);
		}
		if (sd.getDihAngles().size() != d.length) {
			System.out.print("DihSzög lista hossza rossz: desc="
					+ sd.getDihAngles().size() + " d=" + d.length);
			System.exit(1);
		}
		for (int i = 0; i < sd.getDihAngles().size(); i++) {
			System.out.println("moleculeName + próba indul, dih:" + i);
			if (doMapping(i, sd, r, a, d)) {
				break;
			}
		}
	}

	boolean doMapping(int i, StructDescriptor sd,  double[] r, double[] a, double[] d) {
		SphericalToCartesian stc = new SphericalToCartesian();
		ArrayList<Integer> placedAtom = new ArrayList<Integer>();
		int aNum1, aNum2, aNum3, aNum4;
		Double dist;
		// nézzük az elsõ négy atomot az i-edik dihedrál második atomjától
		// kezdve
		// Az elsõ atom marad a 0,0,0-ban
		Vector3D posXYZ = new Vector3D(0, 0, 0);
		StructDescriptor.DihAngleDesc dd = sd.getDihAngles().get(i);
		aNum1 = dd.getA1();
		aNum2 = dd.getA2();
		aNum3 = dd.getA3();
		aNum4 = dd.getA4();
		Atom A1,A2,A3,A4;
		A2=getAtom(aNum2);
		A2.set(posXYZ);
		placedAtom.add(aNum2);
		// A második az X tengely mentén van odább
		A1=getAtom(aNum1);
		dist = findBond(r, sd, aNum1, aNum2);
		if (dist == null)
			return false;
		posXYZ = new Vector3D(dist, 0, 0);
		A1.set(posXYZ);
		placedAtom.add(aNum1);
		// A harmadik az Z mentén is odább van
		A3=getAtom(aNum3);
		dist = findBond(r, sd, aNum2, aNum3);
		if (dist == null)
			return false;
		Double angle = findAngle(a, sd, aNum1, aNum2, aNum3);
		if (angle == null)
			return false;
		posXYZ.set(stc.transformCopy(new double[] { dist, angle, 0d }));   // TODO check 
		A3.set(posXYZ);
		placedAtom.add(aNum3);
		// A negyedik az Y mentén is odább van
		dist = findBond(r, sd, aNum3, aNum4);
		if (dist == null)
			return false;
		angle = findAngle(a, sd, aNum2, aNum3, aNum4);
		if (angle == null)
			return false;
		posXYZ = ZMatrixToCartesian.transform(A1, A2, A3, dist, angle, d[i]/valt); // TODO check
		A4=getAtom(aNum4);
		A4.set(posXYZ);
		placedAtom.add(aNum4);
		int dihedralTried = -1;
		boolean notFound = false;
		// jöhet az összes többi
		while (placedAtom.size() < getAtoms().size()) {
			Integer dih = findNextDihedral(placedAtom, sd, dihedralTried);
//			System.out.println("Placed so far:"	+ placedAtom.size() + placedAtom.toString());
			if (dih != null) {
				dd = sd.getDihAngles().get(dih*Integer.signum(dih));
				if (dih>0){
				aNum1 = dd.getA1();
				aNum2 = dd.getA2();
				aNum3 = dd.getA3();
				aNum4 = dd.getA4();
				}
				else {
				aNum1 = dd.getA4();
				aNum2 = dd.getA3();
				aNum3 = dd.getA2();
				aNum4 = dd.getA1();
				dih=dih*Integer.signum(dih);
				}
				dihedralTried = dih;
//				System.out.println("trying " + dih + ":" + aNum1 + " "+ aNum2 + " " + aNum3 + " " + aNum4);
				dist = findBond(r, sd, aNum3, aNum4);
				if (dist == null) {
					System.out
							.println("NO BOND FOR:" + aNum3 + " " + aNum4);
					continue;
				}
				angle = findAngle(a, sd, aNum2, aNum3, aNum4);
				if (angle == null) {
					System.out.println("NO angle FOR:" + aNum2 + " "
							+ aNum3 + " " + aNum4);
					continue;
				}
				posXYZ = ZMatrixToCartesian.transform(getAtom(aNum1),
						getAtom(aNum2), getAtom(aNum3), dist, angle,
						d[dih]/valt);
				getAtom(aNum4).set(posXYZ);
				placedAtom.add(aNum4);
				notFound = false;
			} else if (dihedralTried != -1 && !notFound) {
				dihedralTried = -1;
				notFound = true;
			} else {
				break;
			}
		}
		if (placedAtom.size() < getAtoms().size()) {
			System.out.println(moleculeName + " " + i
					+ ". dihedrállal sajnos nem sikerült");
			return false;
		} else {
			System.out.println(moleculeName + " " + i
					+ ". dihedrállal sikerült!");
			return true;
		}
	}

	private Double findBond(double[] r, StructDescriptor sd,  int a1, int a2) {
		for (int i = 0; i < sd.getBonds().size(); i++) {
			int ba1 = sd.getBonds().get(i).getA1();
			int ba2 = sd.getBonds().get(i).getA2();
			if ((a1 == ba1 && a2 == ba2) || (a1 == ba2 && a2 == ba1)) {
				return r[i];
			}
		}
		return null;
	}

	private Double findAngle(double[] a, StructDescriptor sd, int a1, int a2, int a3) {
		for (int i = 0; i < sd.getAngles().size(); i++) {
			int ba1 = sd.getAngles().get(i).getA1();
			int ba2 = sd.getAngles().get(i).getA2();
			int ba3 = sd.getAngles().get(i).getA3();
			if ((a2 == ba2)
					&& ((a1 == ba1 && a3 == ba3) || (a1 == ba3 && a3 == ba1))) {
				return (180-a[i])/valt;
			}
		}
		System.out.println("check fail");
		return null;
	}

	private Integer findNextDihedral(ArrayList<Integer> placedAtom, StructDescriptor sd, 
			int dihedralTried) {
		for (int i = dihedralTried + 1; i < sd.getDihAngles().size(); i++) {
			int ba1 = sd.getDihAngles().get(i).getA1();
			int ba2 = sd.getDihAngles().get(i).getA2();
			int ba3 = sd.getDihAngles().get(i).getA3();
			int ba4 = sd.getDihAngles().get(i).getA4();
			if (placedAtom.contains(ba1) && placedAtom.contains(ba2)
					&& placedAtom.contains(ba3) && !placedAtom.contains(ba4)) {
//				System.out.println("check:" + i + "success");
				return i;
			} else if (placedAtom.contains(ba4) && placedAtom.contains(ba3)
					&& placedAtom.contains(ba2) && !placedAtom.contains(ba1)) {
//				System.out.println("check:" + i + "success");
				return -1 * i;
			}
		}
		return null;
	}
}

/*
 * HETATM section of pdb file
 * 
 * HETATM 21 C 0 2.501 -0.658 -2.354 C
 * 
 * COLUMNS DATA TYPE FIELD DEFINITION
 * ----------------------------------------------------------------------- 1 - 6
 * Record name "HETATM" 7 - 11 Integer serial Atom serial number. 13 - 16 Atom
 * name Atom name. 17 Character altLoc Alternate location indicator. 18 - 20
 * Residue name resName Residue name. 22 Character chainID Chain identifier. 23
 * - 26 Integer resSeq Residue sequence number. 27 AChar iCode Code for
 * insertion of residues. 31 - 38 Real(8.3) x Orthogonal coordinates for X. 39 -
 * 46 Real(8.3) y Orthogonal coordinates for Y. 47 - 54 Real(8.3) z Orthogonal
 * coordinates for Z. 55 - 60 Real(6.2) occupancy Occupancy. 61 - 66 Real(6.2)
 * tempFactor Temperature factor. 77 - 78 LString(2) element Element symbol;
 * right-justified. 79 - 80 LString(2) charge Charge on the atom.
 */

/*
 * CONECT section
 * 
 * CONECT 21 20 22 23 24
 * 
 * COLUMNS DATA TYPE FIELD DEFINITION
 * ------------------------------------------------------- 1 - 6 Record name
 * "CONECT" 7 - 11 Integer serial Atom serial number 12 - 16 Integer serial
 * Serial number of bonded atom 17 - 21 Integer serial Serial number of bonded
 * atom 22 - 26 Integer serial Serial number of bonded atom 27 - 31 Integer
 * serial Serial number of bonded atom
 */

