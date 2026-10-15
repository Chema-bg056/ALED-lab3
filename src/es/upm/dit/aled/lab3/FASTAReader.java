package es.upm.dit.aled.lab3;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Reads a FASTA file containing genetic information and allows for the search
 * of specific patterns within these data. The information is stored as an array
 * of bytes that contain nucleotides in the FASTA format. Since this array is
 * usually created before knowing how many characters in the origin FASTA file
 * are valid, an int indicating how many bytes of the array are valid is also
 * stored. All valid characters will be at the beginning of the array.
 * 
 * @author mmiguel, rgarciacarmona
 *
 */
public class FASTAReader {

	protected byte[] content;
	protected int validBytes;

	/**
	 * Creates a new FASTAReader from a FASTA file.
	 * 
	 * @param fileName The name of the FASTA file.
	 */
	public FASTAReader(String fileName) {
		try {
			this.readFile(fileName);
		} catch (IOException e) {
			System.out.println(e.getMessage());
			return;
		}
	}

	/*
	 * Helper method to read from a file. It populates the data array with upper
	 * case version of all the nucleotids found in the file. Throws an IOException
	 * if there is a problem accessing the file or the file is to big to fit in an
	 * array.
	 */
	private void readFile(String fileName) throws IOException {
		File f = new File(fileName);
		FileInputStream fis = new FileInputStream(f);
		DataInput fid = new DataInputStream(fis);
		long len = (int) fis.getChannel().size();
		if (len > Integer.MAX_VALUE) {
			fis.close();
			throw new IOException("The file " + fileName + " is too big. Can't be contained in an array.");
		}
		byte[] content = new byte[(int) len];
		int bytesRead = 0;
		int numRead = 0;
		String line;
		while ((line = fid.readLine()) != null) {
			// Put every character in upper case
			line = line.toUpperCase();
			numRead = line.length();
			byte[] newData = line.getBytes();
			for (int i = 0; i < numRead; i++)
				content[bytesRead + i] = newData[i];
			bytesRead += numRead;
		}
		fis.close();
		this.content = content;
		this.validBytes = bytesRead;
	}

	/**
	 * Provides the data array that contains the complete sequence of nucleotids
	 * extracted from the FASTA file.
	 * 
	 * @return The data array with each nucleotid taking one byte.
	 */
	public byte[] getContent() {
		return content;
	}

	/**
	 * Provides the amount of bytes in the data array that are valid. Since this
	 * array is created before the amount of bytes in the FASTA file that contain
	 * actual nucleotids are know, a worst-case scenario is assumed. So, only
	 * positions between content[0] and content[validBytes -1] have actual genomic
	 * data.
	 * 
	 * @return The number of valid bytes.
	 */
	public int getValidBytes() {
		return validBytes;
	}

	/**
	 * Returns the sequence of nucleotides of the provided size found at the
	 * provided position of the data array. If the initialPos + size is after the
	 * valid bytes of the array, it returns null.
	 * 
	 * @param initialPos The first character of the sequence.
	 * @param size       The length of the sequence.
	 * @return An String representing the sequence.
	 */
	public String getSequence(int initialPos, int size) {
		if (initialPos + size >= validBytes)
			return null;
		return new String(content, initialPos, size);
	}

	/*
	 * Helper method that checks if a pattern appears at a specific position in the
	 * data array. It checks every byte of the pattern one by one. If the pattern
	 * length would need to access a position after the valid bytes of the array, it
	 * throws a new FASTAException to indicate this fact.
	 * 
	 * Returns true if the pattern is present at the given position; false
	 * otherwise.
	 */
	private boolean compare(byte[] pattern, int position) throws FASTAException {
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		boolean match = true;
		for (int i = 0; i < pattern.length; i++) {
			if (pattern[i] != content[position + i]) {
				match = false;
			}
		}
		return match;
	}

	/*
	 * Improved version of the compare method that stops checking elements of the
	 * pattern when one has been found to be different.
	 */
	private boolean compareImproved(byte[] pattern, int position) throws FASTAException {
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		boolean match = true;
		for (int i = 0; i < pattern.length; i++) {
			if (pattern[i] != content[position + i]) {
				return false;
			}
		}
		return match;
		
	}

	/*
	 * Improved version of the compare method that returns the number of bytes in
	 * the pattern that are different from the ones present in the data array at the
	 * given position.
	 * 
	 * Returns the number of characters in the pattern that are different from the
	 * ones present in the indicated position.
	 */
	private int compareNumErrors(byte[] pattern, int position) throws FASTAException {
		int numErrors = 0;
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		for (int i = 0; i < pattern.length; i++) {
			if (pattern[i] != content[position + i]) {
				numErrors++;
			}
		}
		return numErrors;
	}

	/**
	 * Implements a linear search to look for the provided pattern in the data
	 * array. Returns a List of Integers that point to the initial positions of all
	 * the occurrences of the pattern in the data.
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern in the data.
	 */
	public List<Integer> search(byte[] pattern) {
		List<Integer> lista = new ArrayList<Integer>();
			
			for( int i=0; i< content.length; i++) {
				try {
					if (compareImproved(pattern, i)) {
				lista.add(i);
				}
			}
				catch (FASTAException e){
				System.out.println(e.getMessage());
			
		}}
		 return lista;

	}

	/**
	 * Implements a linear search to look for the provided pattern in the data array
	 * but allowing a SNV (Single Nucleotide Variant). In SNV, one nucleotide is
	 * allowed to be different from the pattern. Therefore, this method returns a
	 * List of Integers that point to the initial positions of all the occurrences
	 * of the pattern in the data and all the occurrences of the pattern with one
	 * error in the data
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern (with up to 1 errors) in the data.
	 */
	public List<Integer> searchSNV(byte[] pattern) {
		List<Integer> lista = new ArrayList<Integer>();
		
		for( int i=0; i< content.length; i++) {
			try {
				if (compareNumErrors(pattern, i) <= 1) {
			lista.add(i);
			}
		}
			catch (FASTAException e){
			System.out.println(e.getMessage());
		
	}}
	 return lista;
		
	}

	public static void main(String[] args) {
		long t1 = System.nanoTime();
		FASTAReader reader = new FASTAReader(args[0]);
		if (args.length == 1)
			return;
		System.out.println("Tiempo de apertura de fichero: " + (System.nanoTime() - t1));
		long t2 = System.nanoTime();
		List<Integer> posiciones = reader.searchSNV(args[1].getBytes());
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t2));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado : " + args[1] + " en ningun sitio");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}
