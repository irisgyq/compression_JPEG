function[header,data,len_DC,len_AC,new_G]=comJPEG(G,scalar)

[row,col]=size(G);

%mean-normalization
G=G-128;

%dct-transform each block
dct_G=zeros(size(G));
i=1;
j=1;
while i<=row
    while j<=col
        dct_G(i:i+7,j:j+7)=dct2(G(i:i+7,j:j+7));
        j=j+8;
    end
    i=i+8;
    j=1;
end

%Quantization

%quantization matrix
QM=[16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 
    14 13 16 24 40 57 69 56; 14 17 22 29 51 87 80 62;
    18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 
    49 64 78 87 103 121 120 101; 72 92 95 98 112 100 103 99];
new_QM=scalar*QM;

quan_G=zeros(size(G));

i=1;
j=1;
while i<=row
    while j<=col
        quan_G(i:i+7,j:j+7)=round(dct_G(i:i+7,j:j+7)./new_QM);
        j=j+8;
    end
    i=i+8;
    j=1;
end

%fetch DC terms and AC terms
DC_terms=zeros(1,row*col/64);
AC_terms=zeros(row*col/64,63);

i=1;
j=1;
k=1;
while i<=row
    while j<=col
        DC_terms(1,k)=quan_G(i,j);
        AC_terms(k,1:63)=ZigZagScan(quan_G(i:i+7,j:j+7));
        k=k+1;
        j=j+8;
    end
    i=i+8;
    j=1;
end

bit_AC_terms=zeros(1,row*col*63/64);
k=1;
i=1;
while i<=row*col/64
    bit_AC_terms(1,k:k+62)=AC_terms(i,1:63);
    i=i+1;
    k=k+63;
end

%---------------------------Entropy-coding of DC
%coefficients-----------------------------------
[~,col_DC]=size(DC_terms);

%Calculate DC_residuals
DC_residuals=DC_terms;
for i=2:col_DC
    DC_residuals(1,i)=DC_terms(1,i)-DC_terms(1,i-1);
end
   
%Divide DC_terms into 12 categories
DC_k=zeros(size(DC_residuals));

for i=1:col_DC
    if DC_residuals(1,i)==0
           DC_k(1,i)=0;
    else
        DC_k(1,i)=double(floor(log2(abs(DC_residuals(1,i))))+1);
    end
end

%record sign of the residuals
DC_s=zeros(size(DC_residuals));
for i=1:col_DC
    if DC_residuals(1,i)>0
        DC_s(1,i)=1;
    end
end

%use k-1 bits to represent t
len_t_DC=0;
DC_t=zeros(size(DC_residuals));
for i=1:col_DC
   tmp=DC_k(1,i)-1;
   if tmp>0
           binary=dec2bin(abs(abs(DC_residuals(1,i))-power(2,tmp)),tmp);
           len_t_DC=len_t_DC+length(binary);
           DC_t(1,i)=str2double(binary);
   end
end

%Develop a huffman code for the 12 categories
category_DC=[0,1,2,3,4,5,6,7,8,9,10,11];
p_DC=zeros(size(category_DC));

count=0;
for i=1:12
    for j=1:col_DC
        if double(category_DC(1,i))==DC_k(1,j)
            count=count+1;
        end
    end
    p_DC(i)=count/col_DC;
    count=0;
end

dict_DC=huffmandict(category_DC,p_DC);
DC_code=[];

len_s_DC=0;
len_h_DC=0;
c=cell2mat(dict_DC(:,1));
Stmp="";
max_codeword_DC=0;

for j=1:12
            [~,col_dict_DC]=size(dict_DC{j,2});
            if col_dict_DC>max_codeword_DC
                max_codeword_DC=col_dict_DC;
            end
end
max_codeword_DC=ceil(log2(max_codeword_DC));

l_codeword_DC=zeros(1,12);
len_h0_DC=0;
for j=1:12
            [~,col_dict_DC]=size(dict_DC{j,2});
            len_h0_DC=len_h0_DC+col_dict_DC;
            l_codeword_DC(1,j)=str2double(dec2bin(col_dict_DC,max_codeword_DC));
end

for i=1:col_DC
    for j=1:12
        if DC_k(1,i)==c(j,:)
               [~,col_dict_DC]=size(dict_DC{j,2});
               for k=1:col_dict_DC
                    Stmp=strcat(Stmp,num2str(dict_DC{j,2}(k))); 
                    DC_code(end+1)=str2double(Stmp);
               end
               len_h_DC=len_h_DC+col_dict_DC;
        end
        Stmp="";
    end
    DC_code(end+1)=DC_s(1,i);
    len_s_DC=len_s_DC+1;
    if DC_t(1,i)>0
        DC_code(end+1)=DC_t(1,i);
    end
end

%calculate symbol part of DC_header
symbol_num_DC=ceil(log2(12));
DC_symbol=DC_terms;
for i=1:12
    DC_symbol(1,i)=str2double(dec2bin(i,symbol_num_DC));
end

Stmp="";
%DC_header
DC_header=[];
DC_header(end+8)=str2double(dec2bin(max_codeword_DC,8));
for i=1:12
   DC_header(end+1)=DC_symbol(1,i);
    DC_header(end+1)=l_codeword_DC(1,i);
        if DC_k(1,i)==c(i,:)
            [~,col_dict_DC]=size(dict_DC{i,2});
               for k=1:col_dict_DC
                    Stmp=strcat(Stmp,num2str(dict_DC{i,2}(k))); 
                    DC_header(end+1)=str2double(Stmp);
               end
        end
        Stmp="";
end

len_DC_header=double(8+len_h0_DC+symbol_num_DC*12+12*max_codeword_DC);

len_DC_data=double(len_h_DC+len_s_DC+len_t_DC);

len_DC=len_DC_header+len_DC_data;

%----------------------------------Entropy-coding of AC
%coefficients----------------------
AC_d=[];
x=[];
[~,col_AC]=size(bit_AC_terms);

pos=0;
for i=1:col_AC
    if bit_AC_terms(1,i)~=0
         AC_d(end+1)=i-pos-1;
         x(end+1)=bit_AC_terms(1,i);
         pos=i;
    end
end

[~,col_non_zero_AC]=size(x);

%Divide AC_terms into 10 categories
AC_k=zeros(size(x));

for i=1:col_non_zero_AC
    AC_k(1,i)=double(floor(log2(abs(x(1,i))))+1); 
end

%record sign of the AC_terms
len_s_AC=0;
AC_s=zeros(size(x));
for i=1:col_non_zero_AC
    if x(1,i)>0
        AC_s(1,i)=1;
    end
    len_s_AC=len_s_AC+1;
end

%use k-1 bits to represent t
len_t_AC=0;
AC_t=zeros(size(x));
for i=1:col_non_zero_AC
   tmp=AC_k(1,i)-1;
   if tmp>0
           binary=dec2bin(abs(abs(x(1,i))-power(2,tmp)),tmp);
           len_t_AC=len_t_AC+length(binary);
           AC_t(1,i)=str2double(binary);
   end
end

%calculate r for d
AC_p=zeros(size(x));
AC_r=zeros(size(x));
AC_r_k=cell(size(x));
for i=1:col_non_zero_AC
    AC_p(1,i)=floor(AC_d(1,i)/15);
    AC_r(1,i)=str2double(dec2bin(AC_d(1,i)-15*AC_p(1,i),4));
    AC_r_k(1,i)=cellstr(strcat(dec2bin(AC_d(1,i)-15*AC_p(1,i),4),dec2bin(AC_k(1,i),4)));
    AC_k(1,i)=str2double(dec2bin(AC_k(1,i),4));
end

%calculate the possibility of 152 symbols
symbols=cell(1,152);
k=1;
for i=0:14
    for j=1:10
        s=strcat(dec2bin(i,4),dec2bin(j,4));
        symbols(1,k)=cellstr(s); 
        k=k+1;
    end
end
symbols(151)=cellstr("11110000");
symbols(152)=cellstr("11111111"); %I use "11111111" to represent eob temporarily

p_AC=zeros(1,152);
p_sum=0;
r_count=0;

for i=1:col_non_zero_AC
    p_sum=p_sum+AC_p(1,i);
end  

sum=p_sum+col_non_zero_AC+79*106;


for j=1:150
    for i=1:col_non_zero_AC
        if strcmp(AC_r_k(1,i),symbols(1,j))
            r_count=r_count+1;
        end
    end
    p_AC(j)=r_count/sum;
    r_count=0;
end

p_AC(151)=p_sum/sum;
p_AC(152)=79*106/sum;

dict_AC=huffmandict(symbols,p_AC);
AC_code=[];
len_h_AC=0;
c=cell2mat(dict_AC(:,1));
Stmp="";
max_codeword_AC=0;

for j=1:152
       [~,col_dict_AC]=size(dict_AC{j,2});
       if col_dict_AC>max_codeword_AC
                max_codeword_AC=col_dict_AC;
            
       end
end
max_codeword_AC=ceil(log2(max_codeword_AC));

l_codeword_AC=zeros(1,152);
len_h0_AC=0;
for j=1:152
       
            [~,col_dict_AC]=size(dict_AC{j,2});
            len_h0_AC=len_h0_AC+col_dict_AC;
            l_codeword_AC(1,j)=str2double(dec2bin(col_dict_AC,max_codeword_AC));
       
end

for i=1:col_non_zero_AC
    for j=1:152
        if strcmp(AC_r_k(1,i),c(j,:))
            [~,col_dict_AC]=size(dict_AC{j,2});
            for k=1:col_dict_AC
               Stmp=strcat(Stmp,num2str(dict_AC{j,2}(k))); 
               AC_code(end+1)=str2double(Stmp);
            end
            len_h_AC=len_h_AC+col_dict_AC;
        end
        Stmp="";
    end
    AC_code(end+1)=AC_s(1,i);
    if AC_t(1,i)>0
        AC_code(end+1)=AC_t(1,i);
    end
end

%calculate symbol part of AC_header
symbol_num_AC=8;

Stmp="";
%AC_header
AC_header=[];
AC_header(end+8)=str2double(dec2bin(max_codeword_AC,8));
for i=1:152
   AC_header(end+1)=str2double(symbols(1,i));
    AC_header(end+1)=l_codeword_AC(1,i);
        if strcmp(AC_r_k(1,i),c(j,:))
            [~,col_dict_AC]=size(dict_AC{j,2});
            for k=1:col_dict_AC
               Stmp=strcat(Stmp,num2str(dict_AC{i,2}(k))); 
               AC_header(end+1)=str2double(Stmp);
            end
        end
        Stmp="";
end

len_AC_header=double(8+len_h0_AC+symbol_num_AC*152+152*max_codeword_AC);

len_AC_data=double(len_h_AC+len_s_AC+len_t_AC);

len_AC=len_AC_header+len_AC_data;

data={DC_code,AC_code};
header={DC_header,AC_header};

%Dequantize
dequan_G=zeros(row,col);
i=1;
j=1;
while i<=row
    while j<=col
        dequan_G(i:i+7,j:j+7)=quan_G(i:i+7,j:j+7).*new_QM;
        j=j+8;
    end
    i=i+8;
    j=1;
end

%Apply the inverse DCT transform
idct_G=zeros(row,col);
i=1;
j=1;
while i<=row
    while j<=col
        idct_G(i:i+7,j:j+7)=idct2(dequan_G(i:i+7,j:j+7));
        j=j+8;
    end
    i=i+8;
    j=1;
end

%mean-normalization
new_G=idct_G+128;

    
    


        