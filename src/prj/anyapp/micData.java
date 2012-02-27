package prj.anyapp;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import android.app.Service;
import android.content.Intent;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioTrack;
import android.os.Environment;
import android.os.IBinder;
import android.util.Log;

// window size = 25 ms, shift 10 ms, HMMs and more feature vectors (39 vs 13) need more training data

public class micData extends Service {

	private static boolean isRecording = false; 
	private static final int AUDIO_SAMPLE_FREQ = 16000;

	public Recorder mRecorder;
	public Thread thread;
	
	/** Called when the activity is first created. */
	@Override
	public void onCreate() {
		thread = new Thread	(new Runnable() {
			public void run() {
				record();
			}     
		});
	}

    @Override
	public int onStartCommand(Intent intent, int flags, int startId)
    {
		isRecording = true;
		thread.start();
		return START_STICKY;
    }

    @Override
	public void onDestroy() 
    {
    	isRecording = false;
    	try {
			thread.join();
		} catch (InterruptedException e)  {e.printStackTrace();}
//		play();
    }
    
	
	public void record2(DataOutputStream dos) {
		try {		
			mRecorder = new Recorder();
			mRecorder.StartRecord(dos);
			while(isRecording) {
					if( Thread.interrupted() )
					Log.e("YAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAY", "" + Thread.interrupted());
			}
			mRecorder.StopRecord();		
		} catch (Exception t) {
			Log.e("AudioRecord","Recording Failed - " + t);
		}
		return;
	}
		
	public void record() {
		File file = new File(Environment.getExternalStorageDirectory().getAbsolutePath() + "/MicRecording.pcm");

		// Delete any previous recording.
		if (file.exists())
			file.delete();

		// Create the new file.
		try {
			file.createNewFile();
		} catch (IOException e) {
			throw new IllegalStateException("Failed to create " + file.toString());
		}

		try {
			// Create a DataOuputStream to write the audio data into the saved file.
			OutputStream os = new FileOutputStream(file);
			BufferedOutputStream bos = new BufferedOutputStream(os);
			final DataOutputStream dos = new DataOutputStream(bos);

			// Create a new AudioRecord object to record the audio.
			record2(dos);
			dos.close();

		} catch (Throwable t) {
			Log.e("DataRecord","Recording Failed - " + t);
		}
	}

	public void play() {
		// Get the file we want to playback.
		File file = new File(Environment.getExternalStorageDirectory().getAbsolutePath() + "/MicRecording.pcm");
		// Get the length of the audio stored in the file (16 bit so 2 bytes per short)
		// and create a short array to store the recorded audio.
		int musicLength = (int)(file.length()/2);
		short[] music = new short[musicLength];


		try {
			// Create a DataInputStream to read the audio data back from the saved file.
			InputStream is = new FileInputStream(file);
			BufferedInputStream bis = new BufferedInputStream(is);
			DataInputStream dis = new DataInputStream(bis);

			// Read the file into the music array.
			int i = 0;
			while (dis.available() > 0) {
				music[i] = dis.readShort();
				i++;
			}


			// Close the input streams.
			dis.close();     


			// Create a new AudioTrack object using the same parameters as the AudioRecord
			// object used to create the file.
			AudioTrack audioTrack = new AudioTrack(
											AudioManager.STREAM_MUSIC, 
											AUDIO_SAMPLE_FREQ, 
											AudioFormat.CHANNEL_CONFIGURATION_MONO,
											AudioFormat.ENCODING_PCM_16BIT, 
											musicLength, 
											AudioTrack.MODE_STREAM);
			// Start playback
			audioTrack.play();

			// Write the music buffer to the AudioTrack object
			audioTrack.write(music, 0, musicLength);


		} catch (Throwable t) {
			Log.e("AudioTrack","Playback Failed");
		}
	}

	@Override
	public IBinder onBind(Intent arg0) {
		// TODO Auto-generated method stub
		return null;
	}
}