import java.util.*;
import java.io.*;

public class PeptideGeneration {

	public static void generatePeptides(float expMass, float expPPM, String[] aaName, float[] aaMass) {
		int m = (int)expMass - 18 - 1; //substract mass of H2O and H+
		int[] v = new int[(aaMass.length)];
		for (int i = 0; i < aaMass.length; i++) {
			v[i] = (int)aaMass[i];
		}

		int[][] aaSeqNumb = new int[v.length + 1][m + 1];
		String[][] aaSeqName = new String[v.length + 1][m + 1];
		String[][] aaSeqMass = new String[v.length + 1][m + 1];

		// if m=0 then return empty set to get the m (1 set)
		for (int i = 0; i <= v.length; i++) {
			aaSeqNumb[i][0] = 1;
			aaSeqName[i][0] = "";
			aaSeqMass[i][0] = "";
		}

		// if no aa given, 0 ways to get the m
		for (int j = 1; j <= m; j++) {
			aaSeqNumb[0][j] = 0;
			aaSeqName[0][j] = "";
			aaSeqMass[0][j] = "";
		}

		for (int i = 1; i <= v.length; i++) {
			for (int j = 1; j <= m; j++) {
				if (v[i - 1] <= j) {
					aaSeqNumb[i][j] = aaSeqNumb[i - 1][j] + aaSeqNumb[i][j - v[i - 1]];
					if (aaSeqNumb[i - 1][j] != 0 && aaSeqNumb[i][j - v[i - 1]] != 0) {
						// get Names
						String aa = "";
						String[] tokens = aaSeqName[i][j - v[i-1]].split(";");
						for (int k = 0; k < tokens.length; k++) {
							aa += tokens[k] + aaName[i-1] + ";";
						}
						aaSeqName[i][j] = aaSeqName[i - 1][j] + aa;
						// get Mass
						String mass = "";
						tokens = aaSeqMass[i][j - v[i-1]].split(";");
						for (int k = 0; k < tokens.length; k++) {
							mass += tokens[k] + "-" + v[i-1] + ";";
						}
						aaSeqMass[i][j] = aaSeqMass[i - 1][j] + mass;
					} else if (aaSeqNumb[i][j - v[i - 1]] != 0) {
						// get Names
						String aa = "";
						String[] tokens = aaSeqName[i][j - v[i-1]].split(";");
						for (int k = 0; k < tokens.length; k++) {
							aa += tokens[k] + aaName[i-1] + ";";
						}
						aaSeqName[i][j] = aa;
						// get Mass
						String mass = "";
						tokens = aaSeqMass[i][j - v[i-1]].split(";");
						for (int k = 0; k < tokens.length; k++) {
							if (j - v[i - 1] == 0) {
								mass += v[i-1] + ";";
							} else {
								mass += tokens[k] + "-" + v[i-1] + ";";
							}
						}
						aaSeqMass[i][j] = mass;
					} else if (aaSeqNumb[i - 1][j] != 0) {
						aaSeqName[i][j] = aaSeqName[i - 1][j];
						aaSeqMass[i][j] = aaSeqMass[i - 1][j];
					}
				} else {
					// just copy the value from the top
					aaSeqNumb[i][j] = aaSeqNumb[i - 1][j];
					if (aaSeqNumb[i - 1][j] != 0) {
						aaSeqName[i][j] = aaSeqName[i - 1][j];
						aaSeqMass[i][j] = aaSeqMass[i - 1][j];
					}
				}
			}
		}

		
		// test
		// print out the matrix
		System.out.println();
		System.out.println("# of seqs for each mass:");
		System.out.println("\t\t\tMass");
		System.out.print("\t\t\t");
		for (int j = 0; j <= m; j++) {
			System.out.print(j + "\t");
		}
		System.out.println();
		System.out.print("\t\t");
		for (int j = 0; j <= m; j++) {
			System.out.print("-\t");
		}
		System.out.println();

		System.out.print("Pepties\t");
		System.out.print("0  |\t");
		for (int j = 0; j <= m; j++) {
			System.out.print(aaSeqNumb[0][j] + "\t");
		}
		System.out.println();
		for (int i = 1; i <= v.length; i++) {
			System.out.print("\t\t" + v[i-1] + "  |\t");
			for (int j = 0; j <= m; j++) {
				System.out.print(aaSeqNumb[i][j] + "\t");
			}
			System.out.println();
		}
		

		// print out the peptides
		System.out.println("Predicted Peptides for Mass " + expMass + " with Delta Mass " + expPPM + " PPM (aa mass " + m + "):");
		if (aaSeqNumb[v.length][m] == 0) {
			System.out.println("No solution");
		}
		String[] solPep = aaSeqName[v.length][m].split(";");
		String[] solMass = aaSeqMass[v.length][m].split(";");

		System.out.println();
		System.out.println("8-mer sequences (rounded mass of aa in the sequence):");
		int count = 0;
		for (int k = 0; k < solMass.length; k++) {
			if (solPep[k].length() == 8) {
				System.out.println(count + 1 + " : " + solPep[k] + " (" + solMass[k] + ")");
				count++;
			}
		}			
		System.out.println("There are totally " + count + " predicted 8-mer sequences.");

		System.out.println();
		System.out.println("9-mer sequences (rounded mass of aa in the sequence):");
		count = 0;
		for (int k = 0; k < solMass.length; k++) {
			if (solPep[k].length() == 9) {
				System.out.println(count + 1 + " : " + solPep[k] + " (" + solMass[k] + ")");
				count++;
			}
		}			
		System.out.println("There are totally " + count + " predicted 9-mer sequences.");

		System.out.println();
		System.out.println("10-mer sequences (rounded mass of aa in the sequence):");
		count = 0;
		for (int k = 0; k < solMass.length; k++) {
			if (solPep[k].length() == 10) {
				System.out.println(count + 1 + " : " + solPep[k] + " (" + solMass[k] + ")");
				count++;
			}
		}			
		System.out.println("There are totally " + count + " predicted 10-mer sequences.");

		System.out.println();
		System.out.println("11-mer sequences (rounded mass of aa in the sequence):");
		count = 0;
		for (int k = 0; k < solMass.length; k++) {
			if (solPep[k].length() == 11) {
				System.out.println(count + 1 + " : " + solPep[k] + " (" + solMass[k] + ")");
				count++;
			}
		}			
		System.out.println("There are totally " + count + " predicted 11-mer sequences.");

		System.out.println();
		System.out.println("All sequences (rounded mass of aa in the sequence):");
		count = 0;
		for (int k = 0; k < solMass.length; k++) {
			System.out.println(count + 1 + " : " + solPep[k] + " (" + solMass[k] + ")");
			count++;
		}			
		System.out.println("There are totally " + count + " predicted sequences.");
	}

	public static void main(String[] args) {
		float expMass = Float.parseFloat(args[0]);
		float expPPM = Float.parseFloat(args[1]);
		try {
			Scanner in = new Scanner(new FileReader("aa_mass_test.txt"));
			int i = 0;
			int aaNum = Integer.parseInt(in.nextLine());
			String[] aaName = new String[aaNum];
			float[] aaMass = new float[aaNum];
			while (in.hasNextLine()) {
				String[] tokens = in.nextLine().split("\t");
				aaName[i] = tokens[0];
				aaMass[i] = Float.parseFloat(tokens[1]);
				i++;
			}
			generatePeptides(expMass, expPPM, aaName, aaMass);
		} catch (IOException e) {
     	   e.printStackTrace();
     	}
    }
}