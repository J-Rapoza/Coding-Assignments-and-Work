%Final Project Code
clear
close all

% reading in an images
path= 'C:\Users\bocaj\Documents\WASHU Spring 2023\Numerical Methods\Final Project\data_1\';
img_file= dir(fullfile(path,'*.jpg'));

%create array for image data
img_data=[];

%use a for loop to loop to read through all the images and convert them to grayscale and to colum
%vectors (note that you must change the class of img from uint8 to a double
%in order to do the mean centering part)
for i=1:length(img_file)
    img=imread(fullfile(path,img_file(i).name));
    img=imresize(img,.25);
    img=im2gray(img);
    img=im2double(img);
    img_data=[img_data, img(:)];
end


%mean centering the data
mean_img= mean(img_data, 2); %get the mean of the img_data
mean_center_data = minus(img_data,mean_img); %calculate difference from the mean 

%Apply PCA using SVD
[U,S,V]=svd(mean_center_data, "econ");
V=V'; %transposing the V values
%U is the left singular vector (eigenfaces)
%S is the singular values that tell how important each egienface is in
%reconstructing the faces in "face space"
%V is the right singular vector that gives the linear coefficents needed to
%reconstruct the faces

k=54; %the number of principal components we want to use
A=U(:,1:k)*S(1:k,1:k)*V(1:k,:); %putting back together the principal
%components using only k components


A=plus(A,mean_img); %Adding back the mean 

[height, width]=size(img); %use the size(img) function to figure out the height and width of images

reshape_face= reshape(A,height,width,[]); %reshaping the colum vectors back into something we can visiaulize
implay(reshape_face); %creating a slideshow of the reshaped faces


% Part 2: Analysis of a single image and whether or not it is likely a
% human face
path_2 = 'C:\Users\bocaj\Documents\WASHU Spring 2023\Numerical Methods\Final Project\data_3\img013.jpg'
new_img = imread(path_2); % replace with the path to your new image
new_img = imresize(new_img,.25);
new_img = im2gray(new_img);
new_img = im2double(new_img);
new_img_data=new_img(:);
new_img_data=new_img_data-mean_img;% completing the same steps from aim 1 code
V=V';%transposing V back to its original

projection= U(:,1:k)'*new_img_data*V(1,1:k); %calculating S using the new images column vector data

S_2=diag(projection);% getting the diagonal new S values

S_1=diag(S);%getting the diagonal original S values
threshold =.046; %user defined threshold found by trial and error

if max(abs(S_2)) < max(S_1)*threshold %code to decide and display if the face is human or not
    disp("This is a human face.");
else
    disp("This is most likely not a human face");
end
imshow(path_2);

