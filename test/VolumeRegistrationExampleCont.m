% [214, 264, 13] from moving image should match [99 133 13] in fixed image
% or [129 253 13] --> [50 128 12]
% seems to roughly work!!!

[w.y, w.x, w.z] = Rmoving.intrinsicToWorld(253, 129, 13);
[t.y, t.x, t.z] = transformPointsForward(geomtform, w.y, w.x, w.z);
[f.y, f.x, f.z] = Rfixed.worldToIntrinsic(t.y, t.x, t.z);
f