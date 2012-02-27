package prj.anyapp;

import android.app.Service;
import android.content.Context;
import android.content.Intent;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.IBinder;
import android.util.Log;

public class AccelerometerReader extends Service implements SensorEventListener {
	private float k = (float) 0.99;
	private float accel[] = new float[3];
	private int i = 0;
	private float x,y,z;
	private float sum_x,sum_y,sum_z;
	private SensorManager mySensorManager; // used to acquire sensor access
    private Sensor accSensor; // accelerometer sensor object
    // ...

    @Override
    public void onCreate() {
        // Set SensorManager and acquire a reference to accelerometer sensor
        mySensorManager = (SensorManager)getSystemService(Context.SENSOR_SERVICE);
        accSensor = mySensorManager.getSensorList(Sensor.TYPE_ACCELEROMETER).get(0);
        
        // it may be a good idea to first check the returned Sensor list
        // size to be sure we actually did get an accelerometer reference
        // ...
    }

    // register current SensorEventListener
    @Override
	public int onStartCommand(Intent intent, int flags, int startId)
    {
    	for (int i=0; i<3; i++) {
    		AnyApp.accel[i] = 0;
    	}
        mySensorManager.registerListener(this, accSensor, SensorManager.SENSOR_DELAY_FASTEST);
        return START_STICKY;
    }

    // unregister current SensorEventListener
    @Override
	public void onDestroy()
    {
        mySensorManager.unregisterListener(this);
    }

    // receives accuracy changes form sensor
    public void onAccuracyChanged(Sensor sensor, int accuracy) {

    }


    // receives sensor changes
    public void onSensorChanged(SensorEvent event) {
        if (event.sensor.getType() == Sensor.TYPE_ACCELEROMETER) {
//        	i++;
//        	x = event.values[0];
//        	y = event.values[1];
//        	z = event.values[2];
//
//        	sum_x += x;
//        	sum_y += y;
//        	sum_z += z;
//	
//        	if (i == 10) {
//        		AnyApp.accel[0] = sum_x/10;
//        		AnyApp.accel[1] = sum_y/10;
//        		AnyApp.accel[2] = sum_z/10;
//        	}
        	//ramping from apple sdk 
        	accel[0] = k*accel[0] + (1-k)*event.values[0];
        	accel[1] = k*accel[1] + (1-k)*event.values[1];
        	accel[2] = k*accel[2] + (1-k)*event.values[2];
        	AnyApp.accel[0] = event.values[0] - accel[0];
        	AnyApp.accel[1] = event.values[1] - accel[1];
        	AnyApp.accel[2] = event.values[2] - accel[2];

        }
    }

	@Override
	public IBinder onBind(Intent arg0) {
		return null;
	}

    // ...
}