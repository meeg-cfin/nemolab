function thphz = nemo_cart2sph(xyz)
[thphz(:,1),thphz(:,2),thphz(:,3)]=cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
