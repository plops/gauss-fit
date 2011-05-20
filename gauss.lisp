;; 10.1.1.89.1895.pdf 2000Chen Robust Spot Fitting for Genetic Spot Array Images
;; fulltext.pdf Robust DNA microarray image analysis Norbert Braendle1, Horst Bischof2, Hilmar Lapp3
(declaim (optimize (speed 1) (safety 1) (debug 1)))
(defpackage :gauss
  (:use :cl))
(in-package :gauss)

(deftype num ()
  `single-float)

(defun inv (a b c d)
  "Invert a 2x2 matrix."
  (declare (type num a b c d))
  (let* ((det (- (* a d)
		 (* b c)))
	 (s (/ det)))
    (values (* s d) (* -1 s b)
	    (* -1 s c) (* s a))))
#+nil
(inv 2e0 3e0 3e0 4e0) ;; -4 3 3 -2


(defun .- (a b x y)
  "Subtract two vectors."
  (declare (type num a b x y))
  (values
   (- a x)
   (- b y)))

#+nil
(.- 1 0 1 3)

(defun dot (a b x y)
  "Calculate dot product of vector AB and vector XY."
  (declare (type num a b x y))
  (+ (* a x) (* b y)))
#+nil
(multiple-value-call #'dot (.- 1e0 0e0 2e0 3e0) 3e0 4e0)

(defun mul (a b c d x y)
  "Multiply ABCD matrix with vector XY."
  (declare (type num a b c d x y))
  (values
   (+ (* a x) (* b y))
   (+ (* c x) (* d y))))
#+nil
(multiple-value-call #'mul (inv 1s0 0s0 0s0 1s0) 1s0 0s0)


(defun gauss (px py x y sxx sxy syy)
  (declare (type num px py x y sxx sxy syy))
  (multiple-value-bind (mx my) (.- px py x y)
    (exp (* -.5 (multiple-value-call #'dot mx my 
		     (multiple-value-call #'mul 
		       (inv sxx sxy sxy syy) mx my))))))
#+nil
(gauss .1 .2 .3 .3 .2 .2 .4)


(defun rel-error (z zmodel)
  (declare (type (simple-array num 2) z zmodel)
	   (values num &optional))
  (let* ((z1 (sb-ext:array-storage-vector z))
	 (zmodel1 (sb-ext:array-storage-vector zmodel))
	 (n (length z1))
	 (z0 (let ((sum 0s0))
	       (dotimes (i n)
		 (incf sum (aref z1 i)))
	       (/ sum n)))
	 (top (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (expt (- (aref z1 i) 
				     (aref zmodel1 i))
				  2)))
		sum))
	 (bottom (let ((sum 0s0))
		   (dotimes (i n)
		     (incf sum (expt (- (aref z1 i)
					z0)
				     2)))
		   sum)))
    (/ top bottom)))

(defun calc-unit-position (w &key (h w))
  (declare (type fixnum w h))
  (let ((px (make-array (list h w) :element-type 'num))
	(py (make-array (list h w) :element-type 'num)))
    (dotimes (j h)
      (let ((b (* (- j (floor h 2)) (/ 1s0 h))))
       (dotimes (i w)
	 (let ((a (* (- i (floor w 2)) (/ 1s0 w))))
	   (setf (aref px j i) a
		 (aref py j i) b)))))
    (values px py)))
(defun calc-pixel-position (w &key (h w))
  (declare (type fixnum w h))
  (let ((px (make-array (list h w) :element-type 'num))
	(py (make-array (list h w) :element-type 'num)))
    (dotimes (j h)
      (dotimes (i w)
	(setf (aref px j i) (* 1s0 i)
	      (aref py j i) (* 1s0 j))))
    (values px py)))
(defun estimate-gauss (z px py)
  (declare (type (simple-array num 2) z px py)
	   (values num num num num num &optional))
  (let* ((z1 (sb-ext:array-storage-vector z))
	 (px1 (sb-ext:array-storage-vector px))
	 (py1 (sb-ext:array-storage-vector py))
	 (n (length z1))
	 (total (let ((sum 0s0))
		  (dotimes (i n)
		    (incf sum (aref z1 i)))
		  sum))
	 (mu-x (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (aref px1 i))))
		(/ sum total)))
	 (mu-y (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (aref py1 i))))
		(/ sum total)))
	 (sxx (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i)
			       (expt (- (aref px1 i) mu-x) 2s0))))
		(/ sum total)))
	 (syy (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i)
			       (expt (- (aref py1 i) mu-y) 2s0))))
		(/ sum total)))
	 (sxy (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (* (- (aref px1 i) mu-x) 
					      (- (aref py1 i) mu-y)))))
		(/ sum total))))
    (declare (type num mu-x mu-y)
	     (type (simple-array num 1) px1 py1 z1))
    (values mu-x mu-y sxx sxy syy)))

(defun estimate-amplitude (z)
  (declare (type (simple-array num 2) z)
	   (values num &optional))
  (destructuring-bind (h w) (array-dimensions z)
    (let* ((zmodel (multiple-value-call #'create 
		     w h (do-fit z)))
	   (z1 (sb-ext:array-storage-vector z))
	   (zm1 (sb-ext:array-storage-vector zmodel))
	   (n (length z1))
	   (top (let ((sum 0s0))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (aref zm1 i))))
		 sum))
	  (bottom (let ((sum 0s0))
		    (dotimes (i n)
		      (incf sum (expt (aref zm1 i) 2)))
		    sum)))
      (declare (type (simple-array num 2) zmodel))
     (/ top bottom))))

(defun create (w h x y sxx sxy syy &key (position :pixel))
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy)
	   (values (simple-array num 2) &optional))
  (let ((ar (make-array (list h w) :element-type 'num)))
    (multiple-value-bind (px py)(ecase position
				   (:pixel (calc-pixel-position w :h h))
				   (:unit (calc-unit-position w :h h)))
      (dotimes (j h)       
	(dotimes (i w)
	  (setf (aref ar j i)
		(gauss (aref px j i)
		   (aref py j i)
		   x y sxx sxy syy)))))
    ar))

(defun create-default (w &key (h w) (x 0e0) (y 0e0) (sxx .1s0) (sxy 0s0) (syy sxx) (position :pixel))
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy))
  (create w h x y sxx sxy syy :position position))

(defun scale (q &key (s 9))
  (declare (type (simple-array num 2) q))
  (let* ((r (make-array (array-dimensions q)
			:element-type '(unsigned-byte 8)))
	 (r1 (sb-ext:array-storage-vector r))
	 (q1 (sb-ext:array-storage-vector q)))
    (dotimes (i (length q1))
      (setf (aref r1 i) (min 255 (max 0 (floor (* (+ s .5) (aref q1 i)))))))
    r))

(defun do-fit (img &key (position :pixel))
  (declare (type (simple-array num 2) img))
  (destructuring-bind (h w) (array-dimensions img)
    (multiple-value-bind (px py) (ecase position
				   (:pixel (calc-pixel-position w :h h))
				   (:unit (calc-unit-position w :h h)))
      (estimate-gauss img px py))))

(defun run ()
  (let* ((test (create-default 13 :sxx .04 :sxy .03))
	 (fit (multiple-value-call #'create 13 13 (do-fit test))))
    (format t "~a~%" (list :parameters (multiple-value-list (do-fit test))
			   :amplitude (estimate-amplitude test)
			   :relative-error (rel-error test fit)))
    fit))
#+nil
(scale (run))

#+nil
(scale
 (create-default 13 :x 2s0 :y 2s0 :sxx 4s0 :sxy .03
		 ))


(defun box-muller ()
  "Box Muller random noise generator, polar form"
  ;; http://www.taygeta.com/random/gaussian.html
  (let ((x 0s0)
	(y 0s0)
	(w 2s0)) 
    (loop until (< w 1s0) do
	 (setf x (1- (random 2s0))
	       y (1- (random 2s0))
	       w (+ (* x x) (* y y))))
    (setf w (/ (sqrt (* -2 (log w)))
	       w))
    (values (* x w) (* y w))))

(defun normal-random-number (mean standard-deviation)
  (declare (type num mean standard-deviation))
  (+ (* (box-muller) standard-deviation)
     mean))
#+nil
(normal-random-number 0s0 .1s0)

;; http://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
(defun poisson-random-number (lambda)
  "Knuth's algorithm, only useful for small lambda"
  (declare (type num lambda)
	   (values fixnum &optional))
  (let ((l (exp (- lambda)))
	(k 0)
	(p 1s0))
    (declare (type num l p)
	     (type fixnum k))
    (loop while (< l p) do
	 (incf k)
	 (setf p (* p (random 1s0))))
    (1- k)))
#+nil
(defparameter *bag* (let* ((n 100000)
			   (b (make-array n :element-type 'num)))
		      (dotimes (i n) 
			(setf (aref b i) (normal-random-number 19s0 200s0)))
		      b))

(defun histogram (a &key (n nil))
  (let* ((mi (reduce #'min a))
	 (ma (reduce #'max a))
	 (a-fix (make-array (length a) :element-type 'fixnum))
	 (nn (if n
		 n
		 (floor (- ma mi))))
	 (hist (make-array nn :element-type '(integer 0))))
    (dotimes (i (length a))
      (setf (aref a-fix i)
	    (floor (* (1- nn) (- (elt a i) mi))
		   (- ma mi))))
    (dotimes (i (length a))
      (incf (aref hist (aref a-fix i))))
    (values hist mi ma)))
#+nil
(time 
 (histogram *bag*))


(defun print-histogram (hist mi ma)
  (let ((n (length hist)))
   (dotimes (i n)
     (format t "~a~a~a~%" (+ mi (* 1s0 (- ma mi) (/ i n))) #\Tab (aref hist i)))))
#+nil
(with-open-file (*standard-output* "/dev/shm/o.dat" :if-does-not-exist :create
				   :if-exists :supersede :direction :output)
 (multiple-value-call #'print-histogram
   (histogram *bag* :n 1000)))

(defun write-pgm (filename img)
  (declare (simple-string filename)
           ((array (unsigned-byte 8) 2) img)
           (values null &optional))
  (destructuring-bind (h w)
      (array-dimensions img)
    (declare ((integer 0 65535) w h))
    (with-open-file (s filename
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
      (declare (stream s))
      (format s "P5~%~D ~D~%255~%" w h))
    (with-open-file (s filename 
                       :element-type '(unsigned-byte 8)
                       :direction :output
                       :if-exists :append)
      (let ((data-1d (make-array 
                      (* h w)
                      :element-type '(unsigned-byte 8)
                      :displaced-to img)))
        (write-sequence data-1d s)))
    nil))
(write-pgm "/dev/shm/o.pgm"
	   (scale
	    (create-default 130 :sxx .04 :sxy .03)
	    :s 255))

(print-histogram #(1 2 3) 1 3)
(defun psi (a x)
  (declare (type num x a))
  (if (< (abs x) a)
      (* x (expt (- 1 (expt (/ x a) 2)) 2))
      0s0))

(defun w1 (x)
  (declare (type num x))
  (/ (psi 4s0 x)
     x))



(defun estimate-gauss-robust (z px py old-g 
			      old-mu-x old-mu-y old-sxx old-sxy old-syy old-amp)
  (declare (type (simple-array num 2) z px py old-g)
	   (type num old-mu-x old-mu-y old-sxx old-sxy old-syy))
  (let* ((z1 (sb-ext:array-storage-vector z))
	 (old-g1 (sb-ext:array-storage-vector old-g))
	 (px1 (sb-ext:array-storage-vector px))
	 (py1 (sb-ext:array-storage-vector py))
	 (n (length z1))
	 (sigma (let ((q (make-array n :element-type 'num)))
	      (dotimes (i n)
		(setf (aref q i) (abs (/ (- (* old-amp (aref old-g1 i))
					    (aref z1 i))
					 .6745))))))
	 (total (let ((sum 0s0))
		  (dotimes (i n)
		    (incf sum (aref z1 i)))
		  sum))
	 (mu-x (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (aref px1 i))))
		(/ sum total)))
	 (mu-y (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (aref py1 i))))
		(/ sum total)))
	 (sxx (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (expt (- (aref px1 i) mu-x) 2))))
		(/ sum total)))
	 (syy (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (expt (- (aref py1 i) mu-y) 2))))
		(/ sum total)))
	 (sxy (let ((sum 0s0))
		(dotimes (i n)
		  (incf sum (* (aref z1 i) (* (- (aref px1 i) mu-x) 
					      (- (aref py1 i) mu-y)))))
		(/ sum total))))
    (values mu-x mu-y sxx sxy syy)))