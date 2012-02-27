package prj.anyapp;

import java.io.DataOutputStream;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Iterator;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.AudioRecord.OnRecordPositionUpdateListener;
import android.media.MediaRecorder;
import android.util.Log;

public class Recorder { 
	public static final int AUDIO_SAMPLE_FREQ = 16000; 
	public static final int BUFFER_SIZE = 16000*5; // 5 sec, factor of 1024 
	public  static final int SMALL_BUFFER_SIZE = 320; //10ms
	public  static final int WINDOW_SIZE = 800;
	public  static final int NOTIFICATION_PERIOD = 320; //10ms

	short[] audioBuffer = new short[SMALL_BUFFER_SIZE]; 
	short[] MFCCBuffer = new short[1024]; //MFCC2
	//	float[] MFCCBuffer = new float[1024]; //MFCC
	//	double[] MFCCBuffer = new double[WINDOW_SIZE]; //comirva

	short[][] ab = new short[3][SMALL_BUFFER_SIZE]; 
	//	short[] ab2 = new short[SMALL_BUFFER_SIZE]; 
	//	short[] ab3 = new short[SMALL_BUFFER_SIZE]; 

	public AudioRecord recorder; 
	public micData micdata;
	private int i = 0;
	private DataOutputStream mDos;

	float sampleRate = 8000; 
	int windowSize = 1024; 
	int numberCoefficients = 20; 
	boolean useFirstCoefficient = false; 
	double minFreq = 20; 
	double maxFreq = 4000;
	int numberFilters = 24;
	
	MFCC mfcc = new MFCC(sampleRate, windowSize, numberCoefficients, useFirstCoefficient, minFreq, maxFreq, numberFilters);
	
	
	public Recorder() { 
		try 
		{ 
			// init recorder 
			recorder = new AudioRecord(MediaRecorder.AudioSource.MIC, 
					AUDIO_SAMPLE_FREQ, 
					AudioFormat.CHANNEL_CONFIGURATION_MONO, 
					AudioFormat.ENCODING_PCM_16BIT, 
					BUFFER_SIZE); 
		} 
		catch (IllegalArgumentException e) 
		{ 
			e.printStackTrace(); 
		}
		Log.e("CONSTRUCTOR STATE", Integer.toString(recorder.getState()));
	}

	public OnRecordPositionUpdateListener mNotification 
	= new OnRecordPositionUpdateListener() { 
		public void onPeriodicNotification(AudioRecord arg0) 
		{ 	// read PCM buffer here 	
			int gettingback = recorder.read(audioBuffer, 0, SMALL_BUFFER_SIZE); 
//			Log.e("Recorder Reading", i + " - " + gettingback);
			int index = (i++)%3;
			System.arraycopy(audioBuffer, 0, ab[index], 0, gettingback);
			if (i > 2) {
				if (index == 2) {
					System.arraycopy(ab[0], 0, MFCCBuffer, 0, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[1], 0, MFCCBuffer, SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[2], 0, MFCCBuffer, 2*SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE/2);
				}
				if (index == 0) {
					System.arraycopy(ab[1], 0, MFCCBuffer, 0, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[2], 0, MFCCBuffer, SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[0], 0, MFCCBuffer, 2*SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE/2);
				}
				if (index == 1) {
					System.arraycopy(ab[2], 0, MFCCBuffer, 0, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[0], 0, MFCCBuffer, SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE);
					System.arraycopy(ab[1], 0, MFCCBuffer, 2*SMALL_BUFFER_SIZE, SMALL_BUFFER_SIZE/2);
				}
				float[] y = new float[MFCCBuffer.length];
				for (int i = 0; i < y.length; i++) {
					y[i] = (float) MFCCBuffer[i];
				}
				double tstart = System.currentTimeMillis();
				AnyApp.MFCCData = mfcc.processWindow(y,0);
				double tend = System.currentTimeMillis();
				double t = tend - tstart;
				Log.d("TIME TAKEN = ", Double.toString(t));
//				for (int i=0; i<3; i++) {
//					Log.d("MFCC DATA - ", Double.toString(AnyApp.MFCCData[i]) );
//				}
			}
			
		}							

		public void onMarkerReached(AudioRecord arg0)
		{ 
		} 
	};

	public void StartRecord(DataOutputStream dos) throws InterruptedException 
	{ 	
		mDos = dos;
		int i = recorder.setPositionNotificationPeriod(NOTIFICATION_PERIOD); 
		Log.e("NOTIFICATION PERIOD SET = ",  (i==0)?"Sucess":"ERROR" );
		recorder.setRecordPositionUpdateListener(mNotification);
		recorder.startRecording(); 
		int msamples = recorder.read(audioBuffer, 0, SMALL_BUFFER_SIZE);
		Log.e("INIIAL READ, NUMBER OF READ SAMPLES = ", Integer.toString(msamples));

	} 

	public void StopRecord() 
	{ 
		recorder.stop(); 
	} 
	public void ReleaseRecord() 
	{ 
		recorder.release(); 
	}
}

