package ec.pso;

import ec.pso.Particle;

public class GraphParticle extends Particle {
	private String stringRepresentation;

	public void setStringRepresentation(String str) {
		stringRepresentation = str;
	}

	@Override
	public String toString() {
		return stringRepresentation;
	}
}
