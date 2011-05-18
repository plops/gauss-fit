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
   (- b y))))

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


(defun g (px py x y sxx sxy syy)
  (declare (type num px py x y sxx sxy syy))
  (multiple-value-bind (mx my) (.- px py x y)
    (exp (* -.5 (multiple-value-call #'dot mx my 
		     (multiple-value-call #'mul 
		       (inv sxx sxy sxy syy) mx my))))))
#+nil
(g .1 .2 .3 .3 .2 .2 .4)

(defun create (w &key (h w) (x 0e0) (y 0e0) (sxx .1s0) (sxy 0s0) (syy .1s0))
  (declare (type fixnum w h) 
	   (type num x y sxx sxy syy))
  (let ((ar (make-array (list h w) :element-type 'single-float)))
    (dotimes (j h)
      (let ((b (/ (- j (floor h 2)) h)))
       (dotimes (i w)
	 (let ((a (/ (- i (floor w 2)) w)))
	   (setf (aref ar j i)
		 (g a b x y sxx sxy syy))))))
    ar))

(defun scale (q)
  (declare (type (simple-array single-float 2) q))
  (let* ((r (make-array (array-dimensions q)
			:element-type '(unsigned-byte 8)))
	 (r1 (sb-ext:array-storage-vector r))
	 (q1 (sb-ext:array-storage-vector q)))
    (dotimes (i (length q1))
      (setf (aref r1 i) (floor (* 9.5 (aref q1 i)))))
    r))

(scale
 (create 12))