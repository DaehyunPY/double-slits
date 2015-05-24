character(len = 1) function ch1(no)
	integer :: zero, no
	character(len = 1) :: c1

	zero = ichar('0')

	c1 = char(no +zero)
	ch1 = c1

	return 
end function ch1





character(len = 2) function ch2(no)
	integer :: zero, no, w1, w2
	character(len = 1) :: c1, c2

	zero = ichar('0')

	w2 = int(no/10)
	w1 = no -w2*10

	c1 = char(w1 +zero)
	c2 = char(w2 +zero)

	ch2 = c2//c1

	return 
end function ch2





character(len = 3) function ch3(no)
	integer :: zero, no, w1, w2, w3
	character(len = 1) :: c1, c2, c3

	zero = ichar('0')

	w3 = int(no/100)
	w2 = int((no -w3*100)/10)
	w1 = no -w3*100 -w2*10

	c3 = char(w3 +zero)
	c2 = char(w2 +zero)
	c1 = char(w1 +zero)

	ch3 = c3//c2//c1

	return 
end function ch3





character(len = 4) function ch4(no)
	integer :: zero, no, w1, w2, w3, w4
	character(len = 1) :: c1, c2, c3, c4

	zero = ichar('0')

	w4 = int(no/1000)
	w3 = int((no -w3*100)/100)
	w2 = int((no -w3*100 -w2*10)/10)
	w1 = no -w4*1000 -w3*100 -w2*10

	c4 = char(w4 +zero)
	c3 = char(w3 +zero)
	c2 = char(w2 +zero)
	c1 = char(w1 +zero)

	ch4 = c4//c3//c2//c1

	return 
end function ch4