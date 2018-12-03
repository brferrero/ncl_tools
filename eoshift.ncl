;
; EOSHIFT: end-of-shift
; based on: https://gcc.gnu.org/onlinedocs/gcc-4.8.1/gfortran/EOSHIFT.html
;
;*****************************************************************************;

undef ("eoshift")
function eoshift (NDARRAY:numeric, DIM:integer, SHIFT:integer)
;
; Performs an end-off shift on elements of ARRAY along the dimension of DIM.
; Performs shift to the left if SHIFT is positive and to the right if it is negative.
; 
; Example:
;   V = (1, 2, 3, 4, 5, 6).
;   V = eoshift(V, 1, 2) shifts the elements in V to the left by 2 positions, producing
;   V = (3, 4, 5, 6, 1, 2)
;
local WORK, ds, rank, m, n
begin
    WORK = NDARRAY;
    ds = dimsizes(NDARRAY)
    rank = dimsizes(ds)
    
    if rank.eq.4 then
        ;print ("4D")
        ;A(t,k,m,n)
        m = ds(2)-1
        n = ds(3)-1
        if (SHIFT.gt.0) then
            if (DIM.eq.1) then
                WORK(:,:,0:(m-SHIFT),:) = NDARRAY(:,:,SHIFT:m,:)
                WORK(:,:,(m-SHIFT+1):m,:) = NDARRAY(:,:,0:(SHIFT-1),:)
            else
                WORK(:,:,:,0:(n-SHIFT)) = NDARRAY(:,:,:,SHIFT:n)
                WORK(:,:,:,(n-SHIFT+1):n) = NDARRAY(:,:,:,0:(SHIFT-1))
            end if
        else 
            SSHIFT = abs(SHIFT)
            if (DIM.eq.1) then
                WORK(:,:,SSHIFT:m,:) = NDARRAY(:,:,0:(m-SSHIFT),:)
                WORK(:,:,0:(SSHIFT-1),:) = NDARRAY(:,:,(m-SSHIFT+1):m,:)
            else
                WORK(:,:,:,SSHIFT:n) = NDARRAY(:,:,:,0:(n-SSHIFT))
                WORK(:,:,:,0:(SSHIFT-1)) = NDARRAY(:,:,:,(n-SSHIFT+1):n)
            end if
        end if
    elseif rank.eq.3 then
        ;print ("3D")
        ;A(t,m,n)
        m = ds(1)-1
        n = ds(2)-1
        if (SHIFT.gt.0) then
            if (DIM.eq.1) then
                WORK(:,0:(m-SHIFT),:) = NDARRAY(:,SHIFT:m,:)
                WORK(:,(m-SHIFT+1):m,:) = NDARRAY(:,0:(SHIFT-1),:)
            else
                WORK(:,:,0:(n-SHIFT)) = NDARRAY(:,:,SHIFT:n)
                WORK(:,:,(n-SHIFT+1):n) = NDARRAY(:,:,0:(SHIFT-1))
            end if
        else 
            SSHIFT = abs(SHIFT)
            if (DIM.eq.1) then
                WORK(:,SSHIFT:m,:) = NDARRAY(:,0:(m-SSHIFT),:)
                WORK(:,0:(SSHIFT-1),:) = NDARRAY(:,(m-SSHIFT+1):m,:)
            else
                WORK(:,:,SSHIFT:n) = NDARRAY(:,:,0:(n-SSHIFT))
                WORK(:,:,0:(SSHIFT-1)) = NDARRAY(:,:,(n-SSHIFT+1):n)
            end if
        end if
    elseif rank.eq.2 then
        ;print ("2D")
        ;A(m,n)
        m = ds(0)-1
        n = ds(1)-1
        if (SHIFT.gt.0) then
            if (DIM.eq.1) then
                WORK(0:(m-SHIFT),:) = NDARRAY(SHIFT:m,:)
                WORK((m-SHIFT+1):m,:) = NDARRAY(0:(SHIFT-1),:)
            else
                WORK(:,0:(n-SHIFT)) = NDARRAY(:,SHIFT:n)
                WORK(:,(n-SHIFT+1):n) = NDARRAY(:,0:(SHIFT-1))
            end if
        else 
            SSHIFT = abs(SHIFT)
            if (DIM.eq.1) then
                WORK(SSHIFT:m,:) = NDARRAY(0:(m-SSHIFT),:)
                WORK(0:(SSHIFT-1),:) = NDARRAY((m-SSHIFT+1):m,:)
            else
                WORK(:,SSHIFT:n) = NDARRAY(:,0:(n-SSHIFT))
                WORK(:,0:(SSHIFT-1)) = NDARRAY(:,(n-SSHIFT+1):n)
            end if
        end if
    else 
        ;print("1D")
        ;A(m)
        m = ds-1
        if (SHIFT.gt.0) then
            WORK(0:(m-SHIFT)) = NDARRAY(SHIFT:m)
            WORK((m-SHIFT+1):m) = NDARRAY(0:(SHIFT-1))
        else
            SSHIFT = abs(SHIFT)
            WORK(SSHIFT:m) = NDARRAY(0:(m-SSHIFT))
            WORK(0:(SSHIFT-1)) = NDARRAY((m-SSHIFT+1):m)
        end if
    end if
    return WORK
end

;*****************************************************************************;

procedure run_test()
;
; perform some tests
;
begin
    V = new((/4,3/),float)
    V(0,0) = 1.
    V(0,1) = 2.
    V(0,2) = 3.
    V(1,0) = 4.
    V(1,1) = 5.
    V(1,2) = 6
    V(2,0) = 7.
    V(2,1) = 8.
    V(2,2) = 9.
    V(3,0) = 10.
    V(3,1) = 11.
    V(3,2) = 12.

    print("Original: V(m,n)")
    write_matrix(V,"7f7.1",False)
    ;
    print("Negative Shift: Down - eoshift (V,1,-1)" )
    v = eoshift (V,1,-1)
    write_matrix(v,"7f7.1",False)
    ;
    print("Positive Shift: Up ^ -  eoshift (V,1,1)")
    v = eoshift (V,1,1)
    write_matrix(v,"7f7.1",False)
    ;
    print("Positive Shift: <- - eoshift (V,2,1)")
    v = eoshift (V,2,1)
    write_matrix(v,"7f7.1",False)
    ;
    print("Negative Shift: -> -  eoshift (V,2,-1)")
    v = eoshift (V,2,-1)
    write_matrix(v,"7f7.1",False)

    A = (/1,2,3,4,5,6/)
    a = eoshift (A,1,-2)
    print(A)
    print(a)
end

;begin
; run_test()
;end
