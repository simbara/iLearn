package prj.anyapp;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import android.app.Activity;
import android.content.Intent;
import android.os.Bundle;
import android.os.Handler;
import android.util.Log;
import android.view.View;
import android.view.View.OnClickListener;
import android.widget.Button;
import android.widget.TextView;


public class AnyApp extends Activity implements OnClickListener {
	private static final String TAG = "ServicesDemo";
	private TextView mainText;
	private Button buttonStart;
	private Boolean Running = false;
	private Handler mHandler = new Handler();
	public static float[] accel = new float[3];
	public static float[] MFCCData = new float[12];
//	public static float[] MFCCData = new float[12];
	private int readSensorInterval = 50;
	private ReadSensors ReadSensorsTask = new ReadSensors();
	private BufferedWriter mBufferedWriter;
	long fileOpenTime;


	/** Called when the activity is first created. */
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.main);
		mainText = (TextView) findViewById(R.id.main_text);
		buttonStart = (Button) findViewById(R.id.buttonStart);
		buttonStart.setOnClickListener(this);

	}

	public void onClick(View src) {
		if (!Running) {
			Running = true;
			buttonStart.setText("Stop");
			Log.d(TAG, "onClick: starting srvice");
			startService(new Intent(this, AccelerometerReader.class));
			startService(new Intent(this, micData.class));
			writeToFile();
			mHandler.postDelayed(ReadSensorsTask, readSensorInterval*10);
		} else {
			Running = false;
			buttonStart.setText("Record");
			Log.d(TAG, "onClick: stopping srvice");
			mHandler.removeCallbacks(ReadSensorsTask);
			stopService(new Intent(this, AccelerometerReader.class));
			stopService(new Intent(this, micData.class));
			try {
//				mBufferedWriter.write("</r>");
				mBufferedWriter.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	class ReadSensors implements Runnable 
	{
		public void run() {
			String s;
			mHandler.removeCallbacks(this);
			mHandler.postDelayed(this, readSensorInterval);

			if(mBufferedWriter != null)
			{
//				mainText.setText("\nacc x = "+accel[0] + "\nacc y = "+accel[1] + "\nacc z = " + accel[2] 
//				                 + "\nM1 = "+MFCCData[0] + "\nM2 = "+MFCCData[1] + "\nM3 = "+MFCCData[2]);
				s = System.currentTimeMillis() + ", " + MFCCData[0] + ", " + MFCCData[1] + ", " + MFCCData[2] + ", " + MFCCData[3] + ", " + MFCCData[4] + ", " + MFCCData[5] + ", " + MFCCData[6] + ", " + MFCCData[7] + ", " + MFCCData[8] + ", " + MFCCData[9] + ", " + MFCCData[10] + ", " +MFCCData[11] + ", ";
				s += accel[0] + ", " + accel[1] + ", " + accel[2] + "\n";
//				s += "\t<compass x=\"" + comp[0] + "\" y=\"" + comp[1] + "\" z=\"" + comp[2] + "\" />\n";
				try	{
					mBufferedWriter.write(s);
				} catch (Exception e) {
					Log.e("Error: ", e.getMessage());
				}
			}

		}
	}
	
	public void writeToFile() {
		String sFilename, sPath;
		int	c1;
		File output_xml;
		try
		{
			sPath = "/sdcard/AnyAppData/";
			sFilename = "datafile";
			c1 = 1001;
			do {
				output_xml = new File(sPath + sFilename + c1 + ".txt");
				c1++;
			} while (output_xml.exists());
			(new File(sPath)).mkdirs();
			mBufferedWriter = new BufferedWriter(new FileWriter(output_xml));
			fileOpenTime = System.currentTimeMillis();
		} catch (Exception e) {
			Log.e("Error: ", e.getMessage());
		}
	}
}