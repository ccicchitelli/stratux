package main

import (
	"math"
	"time"
)

var sampleFreq float64 = 500.0
var q0, q1, q2, q3 float64 = 1.0, 0.0, 0.0, 0.0 // estimated orientation quaternion elements with initial conditions
var magX, magY, magZ float64                    // magnetometer measurements

var deltat float64 = 0.002 // sampling period in seconds (shown as 2 ms)
var beta float64 = 2
var zeta float64 = 0
var a_x, a_y, a_z float64                                   // accelerometer measurements
var w_x, w_y, w_z float64                                   // gyroscope measurements in rad/s
var m_x, m_y, m_z float64                                   // magnetometer measurements
var SEq_1, SEq_2, SEq_3, SEq_4 float64 = 1.0, 0.0, 0.0, 0.0 // estimated orientation quaternion elements with initial conditions
var b_x, b_z float64 = 1, 0                                 // reference direction of flux in earth frame
var w_bx, w_by, w_bz float64 = 0.0, 0.0, 0.0                // estimate gyroscope biases error

var attitudeX, attitudeY, attitudeZ, heading float64 = 0.0, 0.0, 0.0, 0.0
var headingHistory [500]float64
var attitudeXhistory [10]float64
var attitudeYhistory [10]float64
var attitudeZhistory [10]float64
var initCount = 0

// Calculates the current heading, optionally compensating for the current attitude
func CalculateHeading() {
	magXtemp := magX
	magYtemp := magY
	magZtemp := magZ
	//these equations account for tilt error
	magXcomp := magXtemp*math.Cos(attitudeY) + magZtemp*math.Sin(attitudeY)
	magYcomp := magXtemp*math.Sin(attitudeX)*math.Sin(attitudeY) + magYtemp*math.Cos(attitudeX) - magZtemp*math.Sin(attitudeX)*math.Cos(attitudeY)
	tempHeading := 180 * math.Atan2(magYcomp, magXcomp) / math.Pi

	if tempHeading < 0 {
		tempHeading += 360
	}

	for i := len(headingHistory) - 1; i > 0; i-- {
		headingHistory[i] = headingHistory[i-1]
	}

	headingHistory[0] = tempHeading

	var total float64 = 0
	for _, value := range headingHistory {
		total += value
	}

	heading = total / float64(len(headingHistory))
}

// Calculates the current attitude represented as X (roll), Y (pitch), and Z (yaw) values as Euler angles.
func CalculateCurrentAttitudeXYZ() {
	var q0a, q1a, q2a, q3a float64
	// q0a = q0
	// q1a = q1
	// q2a = q2
	// q3a = q3
	q0a = SEq_1
	q1a = SEq_2
	q2a = SEq_3
	q3a = SEq_4

	for i := len(attitudeXhistory) - 1; i > 0; i-- {
		attitudeXhistory[i] = attitudeXhistory[i-1]
		attitudeYhistory[i] = attitudeYhistory[i-1]
		attitudeZhistory[i] = attitudeZhistory[i-1]
	}

	attitudeXhistory[0] = math.Atan2(q0a*q1a+q2a*q3a, 0.5-q1a*q1a-q2a*q2a) * 180 / math.Pi
	attitudeYhistory[0] = math.Asin(-2.0*(q1a*q3a-q0a*q2a)) * 180 / math.Pi
	attitudeZhistory[0] = math.Atan2(q1a*q2a+q0a*q3a, 0.5-q2a*q2a-q3a*q3a) * 180 / math.Pi

	var total float64 = 0
	for i := len(attitudeXhistory) - 1; i >= 0; i-- {
		total += attitudeXhistory[i]
	}

	attitudeX = total / float64(len(attitudeXhistory))

	total = 0
	for i := len(attitudeYhistory) - 1; i >= 0; i-- {
		total += attitudeYhistory[i]
	}

	attitudeY = total / float64(len(attitudeYhistory))

	total = 0
	for i := len(attitudeZhistory) - 1; i >= 0; i-- {
		total += attitudeZhistory[i]
	}

	attitudeZ = total / float64(len(attitudeZhistory))
}

// Gets the current attitude and heading.
func GetCurrentAHRS() (float64, float64, float64, float64) {
	return attitudeX, attitudeY, attitudeZ, heading
}

// Gets the current attitude represented as X (roll), Y (pitch), and Z (yaw) values as Euler angles.
func GetCurrentAttitudeXYZ() (float64, float64, float64) {
	return attitudeX, attitudeY, attitudeZ
}

// Gets the current attitude in quaternion form, resulting in no computational load.
func GetCurrentAttitudeQ() (float64, float64, float64, float64) {
	return q0, q1, q2, q3
}

// Gyro input values should be in radians/second, not degrees/second.
// gx, gy, gz: gyroscope values
// ax, ay, az: accelerometer values
// mx, my, mz: magnetometer values
func AHRSupdateOld(gx, gy, gz, ax, ay, az, mx, my, mz float64) {
	initCount++
	if initCount > 5000 { // 10 seconds
		beta = 0.05
	}

	var recipNorm float64                  // vector norm
	var qDot1, qDot2, qDot3, qDot4 float64 // quaternion rate from gyroscopes elements
	var s0, s1, s2, s3 float64             // estimated direction of the gyroscope error
	var hx, hy float64
	var _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3 float64

	// Return if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if (mx == 0.0) && (my == 0.0) && (mz == 0.0) {
		return
	}

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5 * (-q1*gx - q2*gy - q3*gz)
	qDot2 = 0.5 * (q0*gx + q2*gz - q3*gy)
	qDot3 = 0.5 * (q0*gy - q1*gz + q3*gx)
	qDot4 = 0.5 * (q0*gz + q1*gy - q2*gx)

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if !((ax == 0.0) && (ay == 0.0) && (az == 0.0)) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax*ax + ay*ay + az*az)
		ax *= recipNorm
		ay *= recipNorm
		az *= recipNorm

		// Normalise magnetometer measurement
		recipNorm = invSqrt(mx*mx + my*my + mz*mz)
		mx *= recipNorm
		my *= recipNorm
		mz *= recipNorm

		// store magnetometer raw values for later heading calculation
		magX = mx
		magY = my
		magZ = mz

		// Auxiliary variables to avoid repeated arithmetic
		_2q0mx = 2.0 * q0 * mx
		_2q0my = 2.0 * q0 * my
		_2q0mz = 2.0 * q0 * mz
		_2q1mx = 2.0 * q1 * mx
		_2q0 = 2.0 * q0
		_2q1 = 2.0 * q1
		_2q2 = 2.0 * q2
		_2q3 = 2.0 * q3
		_2q0q2 = 2.0 * q0 * q2
		_2q2q3 = 2.0 * q2 * q3
		q0q0 = q0 * q0
		q0q1 = q0 * q1
		q0q2 = q0 * q2
		q0q3 = q0 * q3
		q1q1 = q1 * q1
		q1q2 = q1 * q2
		q1q3 = q1 * q3
		q2q2 = q2 * q2
		q2q3 = q2 * q3
		q3q3 = q3 * q3

		// Reference direction of Earth's magnetic field
		hx = mx*q0q0 - _2q0my*q3 + _2q0mz*q2 + mx*q1q1 + _2q1*my*q2 + _2q1*mz*q3 - mx*q2q2 - mx*q3q3
		hy = _2q0mx*q3 + my*q0q0 - _2q0mz*q1 + _2q1mx*q2 - my*q1q1 + my*q2q2 + _2q2*mz*q3 - my*q3q3
		_2bx = math.Sqrt(hx*hx + hy*hy)
		_2bz = -_2q0mx*q2 + _2q0my*q1 + mz*q0q0 + _2q1mx*q3 - mz*q1q1 + _2q2*my*q3 - mz*q2q2 + mz*q3q3
		_4bx = 2.0 * _2bx
		_4bz = 2.0 * _2bz

		// Gradient decent algorithm corrective step
		s0 = -_2q2*(2.0*q1q3-_2q0q2-ax) + _2q1*(2.0*q0q1+_2q2q3-ay) - _2bz*q2*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx) + (-_2bx*q3+_2bz*q1)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my) + _2bx*q2*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz)
		s1 = _2q3*(2.0*q1q3-_2q0q2-ax) + _2q0*(2.0*q0q1+_2q2q3-ay) - 4.0*q1*(1-2.0*q1q1-2.0*q2q2-az) + _2bz*q3*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx) + (_2bx*q2+_2bz*q0)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my) + (_2bx*q3-_4bz*q1)*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz)
		s2 = -_2q0*(2.0*q1q3-_2q0q2-ax) + _2q3*(2.0*q0q1+_2q2q3-ay) - 4.0*q2*(1-2.0*q1q1-2.0*q2q2-az) + (-_4bx*q2-_2bz*q0)*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx) + (_2bx*q1+_2bz*q3)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my) + (_2bx*q0-_4bz*q2)*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz)
		s3 = _2q1*(2.0*q1q3-_2q0q2-ax) + _2q2*(2.0*q0q1+_2q2q3-ay) + (-_4bx*q3+_2bz*q1)*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx) + (-_2bx*q0+_2bz*q2)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my) + _2bx*q1*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz)
		recipNorm = invSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3) // normalise step magnitude
		s0 *= recipNorm
		s1 *= recipNorm
		s2 *= recipNorm
		s3 *= recipNorm

		// Apply feedback step
		qDot1 -= beta * s0
		qDot2 -= beta * s1
		qDot3 -= beta * s2
		qDot4 -= beta * s3
	}

	// Integrate rate of change of quaternion to yield quaternion
	q0 += qDot1 * (1.0 / sampleFreq)
	q1 += qDot2 * (1.0 / sampleFreq)
	q2 += qDot3 * (1.0 / sampleFreq)
	q3 += qDot4 * (1.0 / sampleFreq)

	// Normalise quaternion
	recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3)
	q0 *= recipNorm
	q1 *= recipNorm
	q2 *= recipNorm
	q3 *= recipNorm
}

func AHRSupdate(w_x, w_y, w_z, a_x, a_y, a_z, m_x, m_y, m_z float64) {
	initCount++
	if initCount > 7500 { // 15 seconds
		beta = 0.05
		zeta = 0.009
	}

	// local system variables
	var norm float64                                                                                                                       // vector norm
	var SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4 float64                                                             // quaternion rate from gyroscopes elements
	var f_1, f_2, f_3, f_4, f_5, f_6 float64                                                                                               // objective function elements
	var J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33, J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64 float64 //objective function Jacobian elements
	var SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4 float64                                                                         // estimated direction of the gyroscope error
	var w_err_x, w_err_y, w_err_z float64                                                                                                  // estimated direction of the gyroscope
	var h_x, h_y, h_z float64

	// normalise the accelerometer measurement
	norm = math.Sqrt(a_x*a_x + a_y*a_y + a_z*a_z)
	a_x /= norm
	a_y /= norm
	a_z /= norm

	// normalise the magnetometer measurement
	norm = math.Sqrt(m_x*m_x + m_y*m_y + m_z*m_z)
	m_x /= norm
	m_y /= norm
	m_z /= norm // computed flux in the earth frame

	// store magnetometer raw values for later heading calculation
	magX = m_x
	magY = m_y
	magZ = m_z

	// auxiliary variables to avoid reapeated calcualtions
	var halfSEq_1 float64 = 0.5 * SEq_1
	var halfSEq_2 float64 = 0.5 * SEq_2
	var halfSEq_3 float64 = 0.5 * SEq_3
	var halfSEq_4 float64 = 0.5 * SEq_4
	var twoSEq_1 float64 = 2.0 * SEq_1
	var twoSEq_2 float64 = 2.0 * SEq_2
	var twoSEq_3 float64 = 2.0 * SEq_3
	var twoSEq_4 float64 = 2.0 * SEq_4
	var twob_x float64 = 2.0 * b_x
	var twob_z float64 = 2.0 * b_z
	var twob_xSEq_1 float64 = 2.0 * b_x * SEq_1
	var twob_xSEq_2 float64 = 2.0 * b_x * SEq_2
	var twob_xSEq_3 float64 = 2.0 * b_x * SEq_3
	var twob_xSEq_4 float64 = 2.0 * b_x * SEq_4
	var twob_zSEq_1 float64 = 2.0 * b_z * SEq_1
	var twob_zSEq_2 float64 = 2.0 * b_z * SEq_2
	var twob_zSEq_3 float64 = 2.0 * b_z * SEq_3
	var twob_zSEq_4 float64 = 2.0 * b_z * SEq_4
	var SEq_1SEq_2 float64
	var SEq_1SEq_3 float64 = SEq_1 * SEq_3
	var SEq_1SEq_4 float64
	var SEq_2SEq_3 float64
	var SEq_2SEq_4 float64 = SEq_2 * SEq_4
	var SEq_3SEq_4 float64
	var twom_x float64 = 2.0 * m_x
	var twom_y float64 = 2.0 * m_y
	var twom_z float64 = 2.0 * m_z

	// compute the objective function and Jacobian
	f_1 = twoSEq_2*SEq_4 - twoSEq_1*SEq_3 - a_x
	f_2 = twoSEq_1*SEq_2 + twoSEq_3*SEq_4 - a_y
	f_3 = 1.0 - twoSEq_2*SEq_2 - twoSEq_3*SEq_3 - a_z
	f_4 = twob_x*(0.5-SEq_3*SEq_3-SEq_4*SEq_4) + twob_z*(SEq_2SEq_4-SEq_1SEq_3) - m_x
	f_5 = twob_x*(SEq_2*SEq_3-SEq_1*SEq_4) + twob_z*(SEq_1*SEq_2+SEq_3*SEq_4) - m_y
	f_6 = twob_x*(SEq_1SEq_3+SEq_2SEq_4) + twob_z*(0.5-SEq_2*SEq_2-SEq_3*SEq_3) - m_z
	J_11or24 = twoSEq_3 // J_11 negated in matrix multiplication
	J_12or23 = 2.0 * SEq_4
	J_13or22 = twoSEq_1 // J_12 negated in matrix multiplication
	J_14or21 = twoSEq_2
	J_32 = 2.0 * J_14or21 // negated in matrix multiplication
	J_33 = 2.0 * J_11or24 // negated in matrix multiplication
	J_41 = twob_zSEq_3    // negated in matrix multiplication
	J_42 = twob_zSEq_4
	J_43 = 2.0*twob_xSEq_3 + twob_zSEq_1 // negated in matrix multiplication
	J_44 = 2.0*twob_xSEq_4 - twob_zSEq_2 // negated in matrix multiplication
	J_51 = twob_xSEq_4 - twob_zSEq_2     // negated in matrix multiplication
	J_52 = twob_xSEq_3 + twob_zSEq_1
	J_53 = twob_xSEq_2 + twob_zSEq_4
	J_54 = twob_xSEq_1 - twob_zSEq_3 // negated in matrix multiplication
	J_61 = twob_xSEq_3
	J_62 = twob_xSEq_4 - 2.0*twob_zSEq_2
	J_63 = twob_xSEq_1 - 2.0*twob_zSEq_3
	J_64 = twob_xSEq_2

	// compute the gradient (matrix multiplication)
	SEqHatDot_1 = J_14or21*f_2 - J_11or24*f_1 - J_41*f_4 - J_51*f_5 + J_61*f_6
	SEqHatDot_2 = J_12or23*f_1 + J_13or22*f_2 - J_32*f_3 + J_42*f_4 + J_52*f_5 + J_62*f_6
	SEqHatDot_3 = J_12or23*f_2 - J_33*f_3 - J_13or22*f_1 - J_43*f_4 + J_53*f_5 + J_63*f_6
	SEqHatDot_4 = J_14or21*f_1 + J_11or24*f_2 - J_44*f_4 - J_54*f_5 + J_64*f_6

	// normalise the gradient to estimate direction of the gyroscope error
	norm = math.Sqrt(SEqHatDot_1*SEqHatDot_1 + SEqHatDot_2*SEqHatDot_2 + SEqHatDot_3*SEqHatDot_3 + SEqHatDot_4*SEqHatDot_4)
	SEqHatDot_1 = SEqHatDot_1 / norm
	SEqHatDot_2 = SEqHatDot_2 / norm
	SEqHatDot_3 = SEqHatDot_3 / norm
	SEqHatDot_4 = SEqHatDot_4 / norm

	// compute angular estimated direction of the gyroscope error
	w_err_x = twoSEq_1*SEqHatDot_2 - twoSEq_2*SEqHatDot_1 - twoSEq_3*SEqHatDot_4 + twoSEq_4*SEqHatDot_3
	w_err_y = twoSEq_1*SEqHatDot_3 + twoSEq_2*SEqHatDot_4 - twoSEq_3*SEqHatDot_1 - twoSEq_4*SEqHatDot_2
	w_err_z = twoSEq_1*SEqHatDot_4 - twoSEq_2*SEqHatDot_3 + twoSEq_3*SEqHatDot_2 - twoSEq_4*SEqHatDot_1

	// compute and remove the gyroscope baises
	w_bx += w_err_x * deltat * zeta
	w_by += w_err_y * deltat * zeta
	w_bz += w_err_z * deltat * zeta
	w_x -= w_bx
	w_y -= w_by
	w_z -= w_bz

	// compute the quaternion rate measured by gyroscopes
	SEqDot_omega_1 = -halfSEq_2*w_x - halfSEq_3*w_y - halfSEq_4*w_z
	SEqDot_omega_2 = halfSEq_1*w_x + halfSEq_3*w_z - halfSEq_4*w_y
	SEqDot_omega_3 = halfSEq_1*w_y - halfSEq_2*w_z + halfSEq_4*w_x
	SEqDot_omega_4 = halfSEq_1*w_z + halfSEq_2*w_y - halfSEq_3*w_x

	// compute then integrate the estimated quaternion rate
	SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat
	SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat
	SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat
	SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat

	// normalise quaternion
	norm = math.Sqrt(SEq_1*SEq_1 + SEq_2*SEq_2 + SEq_3*SEq_3 + SEq_4*SEq_4)
	SEq_1 /= norm
	SEq_2 /= norm
	SEq_3 /= norm
	SEq_4 /= norm

	// compute flux in the earth frame
	SEq_1SEq_2 = SEq_1 * SEq_2 // recompute auxiliary variables
	SEq_1SEq_3 = SEq_1 * SEq_3
	SEq_1SEq_4 = SEq_1 * SEq_4
	SEq_3SEq_4 = SEq_3 * SEq_4
	SEq_2SEq_3 = SEq_2 * SEq_3
	SEq_2SEq_4 = SEq_2 * SEq_4
	h_x = twom_x*(0.5-SEq_3*SEq_3-SEq_4*SEq_4) + twom_y*(SEq_2SEq_3-SEq_1SEq_4) + twom_z*(SEq_2SEq_4+SEq_1SEq_3)
	h_y = twom_x*(SEq_2SEq_3+SEq_1SEq_4) + twom_y*(0.5-SEq_2*SEq_2-SEq_4*SEq_4) + twom_z*(SEq_3SEq_4-SEq_1SEq_2)
	h_z = twom_x*(SEq_2SEq_4-SEq_1SEq_3) + twom_y*(SEq_3SEq_4+SEq_1SEq_2) + twom_z*(0.5-SEq_2*SEq_2-SEq_3*SEq_3)

	// normalise the flux vector to have only components in the x and z
	b_x = math.Sqrt((h_x * h_x) + (h_y * h_y))
	b_z = h_z
}

func isAHRSValid() bool {
	return stratuxClock.Since(mySituation.LastAttitudeTime) < 1*time.Second // If attitude information gets to be over 1 second old, declare invalid.
}

func invSqrt(x float64) float64 {
	return 1.0 / math.Sqrt(x)
}
