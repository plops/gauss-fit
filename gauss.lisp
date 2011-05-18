(declaim (optimize (speed 3) (safety 1) (debug 1)))

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

(defun calc-position (w &key (h w))
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

(defun estimate-gauss (z px py)
  (declare (type (simple-array num 2) z px py))
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



(defun create (w h x y sxx sxy syy)
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy))
  (let ((ar (make-array (list h w) :element-type 'num)))
    (multiple-value-bind (px py) (calc-position w :h h)
      (dotimes (j h)       
	(dotimes (i w)
	  (setf (aref ar j i)
		(g (aref px j i)
		   (aref py j i)
		   x y sxx sxy syy)))))
    ar))

(defun create-default (w &key (h w) (x 0e0) (y 0e0) (sxx .1s0) (sxy 0s0) (syy sxx))
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy))
  (create w h x y sxx sxy syy))

(defun scale (q)
  (declare (type (simple-array num 2) q))
  (let* ((r (make-array (array-dimensions q)
			:element-type '(unsigned-byte 8)))
	 (r1 (sb-ext:array-storage-vector r))
	 (q1 (sb-ext:array-storage-vector q)))
    (dotimes (i (length q1))
      (setf (aref r1 i) (floor (* 9.5 (aref q1 i)))))
    r))

(defun do-fit (img)
  (declare (type (simple-array num 2) img))
  (destructuring-bind (h w) (array-dimensions img)
    (multiple-value-bind (px py) (calc-position w :h h)
      (estimate-gauss img px py))))

(defun run ()
  (let* ((test (create-default 13 :sxx .04 :sxy .03))
	 (fit (multiple-value-call #'create 13 13 (do-fit test))))
    (format t "~a~%" (list :parameters (multiple-value-list (do-fit test))
			   :relative-error (rel-error test fit)))
    fit))
#+nil
(scale (run))

#+nil
(scale
 (create 13 :sxx .04 :sxy .03))
