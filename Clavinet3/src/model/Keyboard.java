package model;

import com.jsyn.Synthesizer;
import com.jsyn.unitgen.LineOut;

import synthesizerBuilder.SynthesizerManager;

public class Keyboard {
	
	private final int numberOfNotes = 12;
	
	private Synthesizer synth;
	private Clavinet[] notes;
	private LineOut lineOut = new LineOut();
	
	public Keyboard() {
		notes = SynthesizerManager.initializeUnits(numberOfNotes);
		synth = SynthesizerManager.initialize(synth, notes, numberOfNotes, lineOut);
	}
	
	public Clavinet getClavinet(int index) {
		return notes[index];
	}
	
	public int clavinetSize() {
		return notes.length;
	}
	
	public void start(int durationSeconds) {
		SynthesizerManager.start(durationSeconds, synth, lineOut);
	}
	
	public void stop() {
		SynthesizerManager.stop(synth);
	}

}
