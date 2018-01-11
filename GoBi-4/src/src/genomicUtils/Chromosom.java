package src.genomicUtils;

import java.util.HashMap;

public class Chromosom {

	HashMap<String, Gene> genesPos = new HashMap<String, Gene>();
	HashMap<String, Gene> genesNeg = new HashMap<String, Gene>();
	
	
	public Chromosom (HashMap<String, Gene> genesPos, HashMap<String, Gene> genesNeg) {
		this.genesPos = genesPos;
		this.genesNeg = genesNeg;
	}
	
	public void addNegGene(Gene g) {
		genesNeg.put(g.getID(), g);
	}
	
	public void addPosGene(Gene g) {
		genesPos.put(g.getID(), g);
	}
	
	public HashMap<String, Gene> genesNeg() {
		return genesNeg;
	}
	
	public HashMap<String, Gene> genesPos() {
		return genesPos;
	}
}
