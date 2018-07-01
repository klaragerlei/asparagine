package Asparagine;

import java.util.ArrayList;

public class StructDescriptor {
      
	class BondDesc{
		private int a1, a2;
		
		public int getA1() {
			return a1;
		}

		public int getA2() {
			return a2;
		}

		public BondDesc(int a1, int a2) {
			this.a1 = a1;
			this.a2 = a2;
		}
	}

	class AngleDesc{
		private int a1,a2,a3;

		public int getA1() {
			return a1;
		}

		public int getA2() {
			return a2;
		}

		public int getA3() {
			return a3;
		}

		public AngleDesc(int a1, int a2, int a3) {
			this.a1 = a1;
			this.a2 = a2;
			this.a3 = a3;
		}
		
	}
	
	class DihAngleDesc{
		private int a1,a2,a3,a4;

		public int getA1() {
			return a1;
		}

		public int getA2() {
			return a2;
		}

		public int getA3() {
			return a3;
		}

		public int getA4() {
			return a4;
		}

		public DihAngleDesc(int a1, int a2, int a3, int a4) {
			this.a1 = a1;
			this.a2 = a2;
			this.a3 = a3;
			this.a4 = a4;
		}
	}
    
	private ArrayList<String> atd=new ArrayList<String>();
	private ArrayList<BondDesc> bd=new ArrayList<BondDesc>();
	private ArrayList<AngleDesc> ad=new ArrayList<AngleDesc>();
	private ArrayList<DihAngleDesc> dd=new ArrayList<DihAngleDesc>();

	public void addAtomType(String at){
		atd.add(at);
	}
	
	public void addBond(int a1, int a2){
		bd.add(new BondDesc(a1, a2));
	}
	
	public void addAngle(int a1, int a2, int a3){
		ad.add(new AngleDesc(a1, a2, a3));
	}

	public void addDihAngle(int a1, int a2, int a3, int a4){
		dd.add(new DihAngleDesc(a1, a2, a3, a4));
	}

	public ArrayList<String> getAtomTypes() {
		return atd;
	}

	public ArrayList<BondDesc> getBonds() {
		return bd;
	}

	public ArrayList<AngleDesc> getAngles() {
		return ad;
	}

	public ArrayList<DihAngleDesc> getDihAngles() {
		return dd;
	}
	
	
	
}
