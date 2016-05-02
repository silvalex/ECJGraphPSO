package ec.pso;

import ec.pso.Particle;

public class GraphParticle extends Particle {
	private static final long serialVersionUID = 1L;

	private String stringRepresentation;
	private double availability;
	private double reliability;
	private double time;
	private double cost;

	public void setStringRepresentation(String str) {
		stringRepresentation = str;
	}

	@Override
	public String toString() {
		return stringRepresentation;
	}

	public void setAvailability(double availability) {
		this.availability = availability;
	}

	public void setReliability(double reliability) {
		this.reliability = reliability;
	}

	public void setTime(double time) {
		this.time = time;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public double getAvailability() {
		return availability;
	}

	public double getReliability() {
		return reliability;
	}

	public double getTime() {
		return time;
	}

	public double getCost() {
		return cost;
	}
}
