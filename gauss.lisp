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
  (declare (type num a b c d)
	   (values num num num num &optional))
  (let* ((det (- (* a d)
		 (* b c)))
	 (s (/ det)))
    (values (* s d) (* -1 s b)
	    (* -1 s c) (* s a))))
#+nil
(inv 2e0 3e0 3e0 4e0) ;; -4 3 3 -2


(defun .- (a b x y)
  "Subtract two vectors."
  (declare (type num a b x y)
	   (values num num &optional))
  (values
   (- a x)
   (- b y)))

#+nil
(.- 1 0 1 3)

(defun dot (a b x y)
  "Calculate dot product of vector AB and vector XY."
  (declare (type num a b x y)
	   (values num &optional))
  (+ (* a x) (* b y)))
#+nil
(multiple-value-call #'dot (.- 1e0 0e0 2e0 3e0) 3e0 4e0)

(defun mul (a b c d x y)
  "Multiply ABCD matrix with vector XY."
  (declare (type num a b c d x y)
	   (values num num &optional))
  (values
   (+ (* a x) (* b y))
   (+ (* c x) (* d y))))
#+nil
(multiple-value-call #'mul (inv 1s0 0s0 0s0 1s0) 1s0 0s0)


(defun gauss (px py x y ixx ixy iyy)
  (declare (type num px py x y ixx ixy iyy)
	   (values num &optional))
  (multiple-value-bind (mx my) (.- px py x y)
    (let ((arg (multiple-value-call #'dot mx my 
				    (mul ixx ixy ixy iyy
					 mx my))))
      (if (< arg 0s0)
	  1s0
	  (exp (* -.5 arg))))))
#+nil
(gauss .1 .2 .3 .3 .2 .2 .4)


(defun rel-error (z zmodel)
  (declare (type (simple-array num 2) z zmodel)
	   (values double-float &optional))
  (let* ((z1 (sb-ext:array-storage-vector z))
	 (zmodel1 (sb-ext:array-storage-vector zmodel))
	 (n (length z1))
	 (z0 (let ((sum 0d0))
	       (dotimes (i n)
		 (incf sum (aref z1 i)))
	       (/ sum n)))
	 (top (let ((sum 0d0))
		(dotimes (i n)
		  (incf sum (expt (- (aref z1 i) 
				     (aref zmodel1 i))
				  2)))
		sum))
	 (bottom (let ((sum 0d0))
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
      (multiple-value-bind (ixx ixy iyx iyy) (inv sxx sxy sxy syy)
       (dotimes (j h)       
	 (dotimes (i w)
	   (setf (aref ar j i)
		 (gauss (aref px j i)
			(aref py j i)
			x y ixx ixy iyy))))))
    ar))

(defun create-default (w &key (h w) (x 0e0) (y 0e0) (sxx .1s0) (sxy 0s0) (syy sxx) (position :pixel))
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy))
  (create w h x y sxx sxy syy :position position))

(defun scale (q &key (s 9))
  #+ nil (declare (type (simple-array num 2) q))
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
    (loop until (< w .9999s0) do
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

#+nil
(write-pgm "/dev/shm/o.pgm"
	   (scale
	    (noise
	     (create-default 80 :x 42s0 :y 32s0 :sxx 40s0 :sxy 10s0))
	    :s 180))

#+nil
(defun copy-array (array)
  (let ((dims (array-dimensions array)))
    (adjust-array
     (make-array dims :element-type (array-element-type array) :displaced-to array)
     dims)))

(defun copy-array (a)
  (let* ((dims (array-dimensions a))
	(b (make-array dims :element-type 'num))
	(a1 (sb-ext:array-storage-vector a))
	(b1 (sb-ext:array-storage-vector b)))
    (dotimes (i (length a1))
      (setf (aref b1 i) (aref a1 i)))
    b))

#+nil
(let* ((n 80)
       (z (create-default n :x 42s0 :y 32s0 :sxx 40s0 :sxy 10s0))
       (noised (noise (copy-array z)))
       (fit (multiple-value-call #'create n n (do-fit z)))
       (fitn (multiple-value-call #'create n n (do-fit noised))))
  (write-pgm "/dev/shm/00data.pgm" (scale z :s 200))
  (write-pgm "/dev/shm/10data-noise.pgm" (scale noised :s 200))
  (write-pgm "/dev/shm/01fit-nonoise.pgm"
	     (scale fit :s 200))
  (write-pgm "/dev/shm/11fit-noise.pgm"
	     (scale fitn :s 200))
  (list (rel-error z fit) (multiple-value-list (do-fit z))
	(rel-error z fitn) (multiple-value-list (do-fit noised))))

;; medsel.pdf approximate median selection
(defun triplet-adjust (a i step)
  (declare (type (simple-array num 1) a)
	   (type fixnum i step))
  (let ((j (+ i step))
	(k (+ i (* 2 step))))
   (if (< (aref a i) (aref a j))
       (if (< (aref a k) (aref a i))
	   (rotatef (aref a i) (aref a j))
	   (when (< (aref a k) (aref a j))
	     (rotatef (aref a j) (aref a k))))
       (if (< (aref a i) (aref a k))
	   (rotatef (aref a i) (aref a j))
	   (when (< (aref a j) (aref a k))
	     (rotatef (aref a j) (aref a k)))))))

(defun selection-sort (a left size step)
  (declare (type (simple-array num 1) a)
	   (type fixnum left size step))
  (loop for i from left below (+ left (* step (- size 1))) by step do
       (let ((min i))
	 (loop for j from (+ i step) below (+ left (* step size)) by step do
	      (when (< (aref a j) (aref a min))
		(setf min j)))
	 (rotatef (aref a i) (aref a min)))))

(defun approximate-median (a &key (size (length a)) (threshold 8))
  "THRESHOLD is the maximum size of the array, that is sorted. Bigger
values result in better estimate of the median at the expense of
processing time."
  (declare (type (simple-array num 1) a)
	   (type fixnum size threshold))
  (let ((left-to-right nil)
	(left 0)
	(step 1))
    (loop while (< threshold size) do
	 (setf left-to-right (not left-to-right))
	 (let* ((rem (mod size 3))
	       (i (if left-to-right
		      left
		      (+ left (* (+ 3 rem) step)))))
	   (dotimes (j (1- (floor size 3)))
	     (triplet-adjust a i step)
	     (incf i (* 3 step)))
	   (if left-to-right
	       (incf left step)
	       (setf i left
		     left (+ left (* (1+ rem) step))))
	   (selection-sort a i (+ 3 rem) step)
	   (when (= rem 2)
	     (if left-to-right
		 (rotatef (aref a (+ i step)) (aref a (+ i (* 2 step))))
		 (rotatef (aref a (+ i (* 2 step))) (aref a (+ i (* 3 step))))))
	   (setf step (* 3 step)
		 size (floor size 3))))
    (selection-sort a left size step)
    (aref a (+ left (* step (floor (1- size) 2))))))

(defun median (a)
  (let* ((a1 (sb-ext:array-storage-vector a))
	 (s (sort a1 #'<))
	 (n (length a1)))
    (if (oddp n)
	(elt s (floor n 2))
	(* .5 (+ (elt s (floor n 2))
		 (elt s (1+ (floor n 2))))))))

#+nil
(let* ((a (create-default 400 :position :unit))
       (b (create-default 400 :position :unit))
       (a1 (sb-ext:array-storage-vector a)))
  (list (time (approximate-median a1 :threshold 2000))
	(time (median b))))

(defun noise (a &key (type :gaussian))
  (let ((a1 (sb-ext:array-storage-vector a)))
    (dotimes (i (length a1))
      (setf (aref a1 i) (+ (aref a1 i) (normal-random-number 0s0 .001s0))))
    a))

#+nil
(defun psi (a x)
  (declare (type num x a))
  (if (< (abs x) a)
      (* x (expt (- 1 (expt (/ x a) 2)) 2))
      0s0))

#+nil
(defun w1 (x)
  (declare (type num x))
  (if (= x 0s0)
    0s0
    (/ (psi 4s0 x)
       x)))

(defun w1 (x)
  (declare (type num x))
  (let ((a 4s0))
    (if (<= (abs x) a)
	(expt (- 1 (expt (/ x a) 2)) 2)
	0s0)))


(defun draw-covariance-ellipse (w h px py mu-x mu-y sxx sxy syy)
  (let* ((c (expt 1.6449 2))
	 (q (make-array (list h w) :element-type '(unsigned-byte 8))))
    (multiple-value-bind (ixx ixy iyx iyy) (inv sxx sxy sxy syy)
     (dotimes (j h)
       (dotimes (i w)
	 (multiple-value-bind (vx vy) (.- (aref px j i) (aref py j i)
					  mu-x mu-y)
	   (let* ((r (multiple-value-call #'dot vx vy
					  (mul ixx ixy iyx iyy
					       vx vy))))
	     (when (< r c)
	       (setf (aref q j i) 1)))))))
    q))


#+nil
(progn
 (defparameter *ellipse*
   (let ((h 120)
	 (w 120))
     (multiple-value-call #'draw-covariance-ellipse
       w h
       (calc-pixel-position w :h h) 
       (do-fit
	   (create-default w :h h :x 42s0 :y 32s0 
			   :sxx 40s0 :sxy 10s0 :position :pixel)))))
 (write-pgm "/dev/shm/ellipse.pgm" *ellipse*))

(defun .+ (a b)
  (let* ((c (make-array (array-dimensions a) :element-type 'num))
	 (a1 (sb-ext:array-storage-vector a))
	 (b1 (sb-ext:array-storage-vector b))
	 (c1 (sb-ext:array-storage-vector c)))
    (dotimes (i (length c1))
      (setf (aref c1 i) (+ (aref a1 i) (aref b1 i))))
    c))

(defun median-absolute-deviation (a)
  (declare (type (simple-array num 1) a))
  (let* ((m (median a))
	 (abd (make-array (length a) :element-type 'num)))
    (map-into abd #'(lambda (x) (abs (- x m)))
	      a)
    (median abd)))

(defun estimate-gauss-robust (z px py old-fit 
			      old-mu-x old-mu-y 
			      old-sxx old-sxy old-syy old-amp)
  (declare (type (simple-array num 2) z px py old-fit)
	   (type num old-mu-x old-mu-y old-sxx old-sxy old-syy old-amp))
  (destructuring-bind (h w) (array-dimensions z)
   (let* ((dims (array-dimensions z))
	  (z1 (sb-ext:array-storage-vector z))
	  (old-fit1 (sb-ext:array-storage-vector old-fit))
	  (px1 (sb-ext:array-storage-vector px))
	  (py1 (sb-ext:array-storage-vector py))
	  (n (length z1))
	  (residual (let ((q (make-array dims :element-type 'num)))
		      (dotimes (j h)
			(dotimes (i w)
			  (setf (aref q j i) (- (aref z j i)
						(aref old-fit j i)))))
		      q))
	  (residual-abs (let ((q (make-array dims :element-type 'num)))
			  (dotimes (j h)
			    (dotimes (i w)
			      (setf (aref q j i) (abs (aref residual j i)))))
			  q))
	  (ellipse (draw-covariance-ellipse w h px py old-mu-x old-mu-y
					    old-sxx old-sxy old-syy))
	  (inside-points (reduce #'+ (sb-ext:array-storage-vector ellipse)))
	  (inside-region (let ((q (make-array inside-points :element-type 'num))
			       (k 0))
			   (dotimes (j h)
			    (dotimes (i w)
			      (when (= 1 (aref ellipse j i)) 
				(setf (aref q k) (aref residual j i))
				(incf k))))
			   q))
	  (sigma (/ (median-absolute-deviation inside-region) 
		.6745)
	    #+nil (let ((m (median inside-region))) ;; replace residual with residual-abs
	       (/ m .6745)))
	  (weight (let* ((q (make-array dims :element-type 'num)))
		    (when (<= sigma 0s0)
		      (error "sigma not positive ~A" sigma))
		    (dotimes (j h)
		      (dotimes (i w)
			(setf (aref q j i) (w1 (/ (aref residual j i) sigma)))))
		    q))
	  (weight1 (sb-ext:array-storage-vector weight))
	  (total (let ((sum 0s0))
		   (dotimes (i n)
		     (incf sum (aref z1 i)))
		   sum))
	  (total-weight (let ((sum 0s0))
			  (dotimes (i n)
			    (incf sum (* (aref weight1 i) 
					 (aref z1 i))))
			  sum))
	  (mu-x (let ((sum 0s0))
		 (when (= 0s0 total-weight)
		   (error "total-weight is zero"))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (aref px1 i) (aref weight1 i))))
		 (/ sum total-weight)))
	  (mu-y (let ((sum 0s0))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (aref py1 i) (aref weight1 i))))
		 (/ sum total-weight)))
	  (sxx (let ((sum 0s0))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (expt (* (aref weight1 i) 
						     (- (aref px1 i) mu-x)) 2))))
		 (/ sum total)))
	  (syy (let ((sum 0s0))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (expt (* (aref weight1 i) 
						     (- (aref py1 i) mu-y)) 2))))
		 (/ sum total)))
	  (sxy (let ((sum 0s0))
		 (dotimes (i n)
		   (incf sum (* (aref z1 i) (* (- (aref px1 i) mu-x) 
					       (- (aref py1 i) mu-y)
					       (expt (aref weight1 i) 2)))))
		 (/ sum total)))
	  (a (let ((sum 0s0)
		   (bottom 0s0))
	       (dotimes (i n)
		 (let* ((g (/ (aref old-fit1 i) old-amp))
			(v (* g (aref weight1 i))))
		   (incf bottom (* g v))
		   (incf sum (* (aref z1 i) v))))
	       (/ sum bottom))))
     (defparameter *residual* residual)
     (defparameter *ellipse* ellipse)
     (defparameter *weight* weight)
     (defparameter *sigma* sigma)
     (values mu-x mu-y sxx sxy syy a))))

(defun create-robust (w h x y sxx sxy syy a &key (position :pixel))
  (declare (type fixnum w h) 
	   (type num x y a sxx sxy syy)
	   (values (simple-array num 2) &optional))
  (let ((ar (make-array (list h w) :element-type 'num)))
    (multiple-value-bind (px py)(ecase position
				  (:pixel (calc-pixel-position w :h h))
				  (:unit (calc-unit-position w :h h)))
      (multiple-value-bind (ixx ixy iyx iyy) (inv sxx sxy sxy syy)
       (dotimes (j h)       
	 (dotimes (i w)
	   (setf (aref ar j i)
		 (* a (gauss (aref px j i)
			     (aref py j i)
			     x y ixx ixy iyy)))))))
    ar))

(defun v+ (a b)
  (declare (type (simple-array num 2) a b)
	   (values (simple-array num 2) &optional))
  (let ((c (make-array (array-dimensions a)
		       :element-type 'num)))
    (map-into (sb-ext:array-storage-vector c)
	      #'(lambda (x y) (declare (type num x y)) 
		  (+ x y)) 
	      (sb-ext:array-storage-vector a)
	      (sb-ext:array-storage-vector b))
    c))

(defun v* (s a)
  (declare (type (simple-array num 2) a)
	   (type num s)
	   (values (simple-array num 2) &optional))
  (let ((c (make-array (array-dimensions a)
		       :element-type 'num)))
    (map-into (sb-ext:array-storage-vector c)
	      #'(lambda (x) (declare (type num x)) 
		  (* s x)) 
	      (sb-ext:array-storage-vector a))
    c))

#+nil
(let* ((a (create-default 100 :x 50s0 :y 50s0 :sxx 10s0))
       (b (create-default 100 :x 40s0 :y 50s0 :sxx 40s0))
       (c (v+ a b)))
  (write-pgm "/dev/shm/fasl.pgm"
	     (scale c :s 230)))

(defun write-normalized-pgm (fn a)
  (write-pgm fn
	     (scale a
		    :s (/ 255s0 (reduce #'max (sb-ext:array-storage-vector a))))))

(defun gnuplot (&rest rest)
  (when (listp rest)
   (with-open-file (s "/dev/shm/o.dat" :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
     (let* ((el (elt rest 0))
	    (first-array (if (listp el)
			     (second el)
			     el)))
       (dotimes (i (length first-array))
	(format s "~a ~{~a ~}~%" 
		i 
		(loop for e in rest collect
		     (if (listp e)
			 (destructuring-bind (title array) e
			   (aref array i))
			 (aref e i)))))))
   (with-open-file (s "/dev/shm/o.gp" :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
     (format s "plot ~:{\"/dev/shm/o.dat\" u 1:~A w l title ~s~a ~}~% pause -1"
	     (loop for i below (length rest) collect
		  (let ((e (elt rest i)))
		    (if (listp e)
			(destructuring-bind (title array) e
			  (list (+ 2 i)
				title
				(if (= i (1- (length rest)))  #\; #\,)))
			(list (+ 2 i)
			      (format nil "~a" (+ 2 i))
			      (if (= i (1- (length rest)))  #\; #\,)))))
	     ))))

(defun horline (a)
  (destructuring-bind (h w) (array-dimensions a)
   (let ((r (make-array w :element-type 'num)))
     (dotimes (i w)
       (setf (aref r i) (* 1s0 (aref a (floor h 2) i))))
     r)))

#+nil
(let* ((w 64)
       (h 64)
       (x 20s0) (y 32s0) (sxx 6s0) (sxy 0s0) (syy sxx) (amp 1s0)
       (z (create-default w :h h :x x :y y :sxx sxx :sxy sxy))
       (za (v+ z (create-robust w h (+ x 32s0) y sxx 0s0 syy .25s0)))
       (fit ()))
  (write-pgm "/dev/shm/00z.pgm" (scale z :s 200))
  (write-pgm "/dev/shm/01za.pgm" (scale za :s 180))
  (multiple-value-bind (px py) (calc-pixel-position w :h h)
    (multiple-value-bind (x y sxx sxy syy) (estimate-gauss za px py)
      (let ((amp (estimate-amplitude za)))
	(dotimes (k 4)
	  (setf fit 
		#-nil (multiple-value-call #'create-robust w h (do-fit za) (estimate-amplitude za))
		#+nil (create-robust w h x y sxx sxy syy amp))
	  
	  (multiple-value-setq (x y sxx sxy syy amp)
	    (multiple-value-call #'estimate-gauss-robust
	      za
	      px py
	      fit
	      x y sxx sxy syy amp))
	  (write-normalized-pgm (format nil "/dev/shm/01residual-~3,'0d.pgm" k)
				*residual*)
	  (write-normalized-pgm (format nil "/dev/shm/02ellipse-~3,'0d.pgm" k)
				*ellipse*)
	  (write-normalized-pgm (format nil "/dev/shm/03weight-~3,'0d.pgm" k)
				*weight*)
	  (write-pgm (format nil "/dev/shm/04fit-~3,'0d.pgm" k)
				(scale fit :s 200))
	  (gnuplot `("za" ,(horline za))
		   `("fit" ,(horline fit)) 
		   `("residual" ,(horline *residual*))
		   `("ellipse" ,(horline *ellipse*))
		   `("weight" ,(horline *weight*))
		   `("sigma" ,(make-array w :element-type 'num
				     :initial-contents (loop for i below w collect *sigma*))))
	  (format t "~a~%" (list k (rel-error za fit) x y sxx syy amp)))))))
#+nil
(let* ((w 64)
       (h 64)
       (x 20s0) (y 32s0) (sxx 6s0) (sxy 0s0) (syy sxx) (amp 1s0)
       (z (v+
	   (create-robust w h (+ x 8) y sxx 0s0 syy .2s0)
	   (create-default w :h h :x x :y y :sxx sxx :sxy sxy)))
       (fit (multiple-value-call #'create-robust w h (do-fit z) (estimate-amplitude z))))
  (gnuplot (horline fit) (horline z)))

#+nil
(let* ((n 200)
       (z (create-default n :x 142s0 :y 132s0 :sxx 110s0 :sxy 10s0))
       (za (v+ z
	       (create-default n :x 142s0 :y 102s0 :sxx 89s0 :sxy 10s0)))
       (noised (noise (copy-array z)))
       (fit (multiple-value-call #'create n n (do-fit z)))
       (fitn (multiple-value-call #'create n n (do-fit noised)))
       (fiti (multiple-value-call #'create-robust n n 
				  (multiple-value-call #'estimate-gauss-robust
				    z 
				    (calc-pixel-position n)
				    z 140s0 112s0 100s0 .1s0 100s0
				    (estimate-amplitude z)))))
  (write-pgm "/dev/shm/00data.pgm" (scale z :s 200))
  (write-pgm "/dev/shm/10data-noise.pgm" (scale noised :s 200))
  (write-pgm "/dev/shm/01fit-nonoise.pgm"
	     (scale fit :s 200))
  (write-pgm "/dev/shm/11fit-noise.pgm"
	     (scale fitn :s 200))
  (write-pgm "/dev/shm/20fit-robust.pgm"
 	     (scale fitn :s 200))
  (write-pgm "/dev/shm/30fit-two.pgm"
	     (scale fiti  :s 200))
  (list (rel-error z fit) (multiple-value-list (do-fit z))
	(rel-error z fitn) (multiple-value-list (do-fit noised))
	(rel-error z fiti) #+nil (multiple-value-list
			    (let ((amp (estimate-amplitude noised)))
			      (multiple-value-bind (px py) (calc-pixel-position n)
				(multiple-value-call #'estimate-gauss-robust noised px py 
						     fit (do-fit z) amp))))))