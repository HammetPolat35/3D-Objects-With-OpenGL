#CENG-487 Assignment-4
# 320201105 HAMMET POLAT
#05/2025

from OpenGL.GL import *
from vec3d import Vec3d
from cam3d import Cam3D
from math import pi,sin,cos,sqrt,acos
cam = Cam3D()
class Geo:
    @staticmethod
    def makeTriangle(a):
        from object3d import Object3D
        arrTri=[]
        dis = 1 / 2 ** a
        for j in range(0, 2 ** a):
            for i in range(0, j + 1):
                v1 = [
                    Vec3d([1.0-dis*j*2+i*dis, -1.0+i*dis, 0.0+i*dis, 1.0]),#define pivot point
                    Vec3d([1.0-dis*(j+1)*2+i*dis, -1.0+i*dis, 0.0+i*dis, 1.0]),#increase j by 1 from pivot
                    Vec3d([1.0-dis*(j+1)*2+(i+1)*dis, -1.0+(i+1)*dis, 0.0+(1+i)*dis, 1.0])#increase j and i by 1 from pivot
                ]
                v2 = [
                    Vec3d([0-dis*j+dis*i,1-2*dis*j+dis*i,0+dis*i, 1.0]),
                    Vec3d([0-dis*(j+1)+dis*i,1-2*dis*(j+1)+dis*i,0+dis*i, 1.0]),
                    Vec3d([0-dis*(j+1)+dis*(i+1),1-2*dis*(j+1)+dis*(i+1),0+dis*(i+1), 1.0])
                ]
                v3 = [
                    Vec3d([0+j*dis-i*dis,1-2*dis*j+dis*i,0+dis*i,1]),
                    Vec3d([0+(j+1)*dis-i*dis,1-2*dis*(j+1)+dis*i,0+dis*i,1]),
                    Vec3d([0+(j+1)*dis-(i+1)*dis,1-2*dis*(j+1)+dis*(i+1),0+dis*(i+1),1])
                ]
                v4 = [
                    Vec3d([0+dis*j-2*dis*i,1-2*dis*j,0,1]),
                    Vec3d([0+dis*(j+1)-2*dis*i,1-2*dis*(j+1),0,1]),
                    Vec3d([0+dis*(j+1)-2*dis*(i+1),1-2*dis*(j+1),0,1])
                ]

                arrTri.extend([Object3D(v1),Object3D(v2),Object3D(v3),Object3D(v4)])

        for j in range(0, 2 ** a):
            for i in range(0, j):
                v1 = [
                    Vec3d([1.0-dis*j*2+i*dis, -1.0+i*dis, 0.0+i*dis, 1.0]),#define pivot
                    Vec3d([1.0-dis*j*2+(i+1)*dis, -1.0+(i+1)*dis, 0.0+(i+1)*dis, 1.0]),#increase i
                    Vec3d([1.0-dis*(j+1)*2+(i+1)*dis, -1.0+(i+1)*dis, 0.0+(1+i)*dis, 1.0])#increase i and j
                ]
                v2 = [
                    Vec3d([0-dis*j+dis*i,1-2*dis*j+dis*i,0+dis*i, 1.0]),
                    Vec3d([0-dis*j+dis*(i+1),1-2*dis*j+dis*(i+1),0+dis*(i+1), 1.0]),
                    Vec3d([0-dis*(j+1)+dis*(i+1),1-2*dis*(j+1)+dis*(i+1),0+dis*(i+1), 1.0])
                ]
                v3 = [
                    Vec3d([0+j*dis-i*dis,1-2*dis*j+dis*i,0+dis*i,1]),
                    Vec3d([0+j*dis-(i+1)*dis,1-2*dis*j+dis*(i+1),0+dis*(i+1),1]),
                    Vec3d([0+(j+1)*dis-(i+1)*dis,1-2*dis*(j+1)+dis*(i+1),0+dis*(i+1),1])
                ]
                v4 = [
                    Vec3d([0+dis*j-2*dis*i,1-2*dis*j,0,1]),
                    Vec3d([0+dis*j-2*dis*(i+1),1-2*dis*j,0,1]),
                    Vec3d([0+dis*(j+1)-2*dis*(i+1),1-2*dis*(j+1),0,1])
                ]
                arrTri.extend([Object3D(v1), Object3D(v2), Object3D(v3), Object3D(v4)])

        return arrTri



    @staticmethod
    def makeCube(a):
        from object3d import Object3D
        arrSqu = []

        dis = 2 / (2 ** a)

        for j in range(2 ** a):
            for i in range(2 ** a):
                # Square 1 (top square)
                squVer1 = [
                    Vec3d([-1 + i * dis, 1 - j * dis, 1.0, 1.0]),
                    Vec3d([-1 + i * dis, 1 - (j + 1) * dis, 1.0, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1 - (j + 1) * dis, 1.0, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1 - j * dis, 1.0, 1.0])
                ]

                # Square 2 (bottom square)
                squVer2 = [
                    Vec3d([-1 + i * dis, 1 - j * dis, -1.0, 1.0]),
                    Vec3d([-1 + i * dis, 1 - (j + 1) * dis, -1.0, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1 - (j + 1) * dis, -1.0, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1 - j * dis, -1.0, 1.0])
                ]

                # Square 3 (right square)
                squVer3 = [
                    Vec3d([1, 1 - j * dis, -1 + i * dis, 1.0]),
                    Vec3d([1, 1 - (j + 1) * dis, -1 + i * dis, 1.0]),
                    Vec3d([1, 1 - (j + 1) * dis, -1 + (i + 1) * dis, 1.0]),
                    Vec3d([1, 1 - j * dis, -1 + (i + 1) * dis, 1.0])
                ]

                # Square 4 (left square)
                squVer4 = [
                    Vec3d([-1, 1 - j * dis, -1 + i * dis, 1.0]),
                    Vec3d([-1, 1 - (j + 1) * dis, -1 + i * dis, 1.0]),
                    Vec3d([-1, 1 - (j + 1) * dis, -1 + (i + 1) * dis, 1.0]),
                    Vec3d([-1, 1 - j * dis, -1 + (i + 1) * dis, 1.0])
                ]

                # Square 5 (front square)
                squVer5 = [
                    Vec3d([-1 + i * dis, 1, 1 - j * dis, 1.0]),
                    Vec3d([-1 + i * dis, 1, 1 - (j + 1) * dis, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1, 1 - (j + 1) * dis, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, 1, 1 - j * dis, 1.0])
                ]

                # Square 6 (back square)
                squVer6 = [
                    Vec3d([-1 + i * dis, -1, 1 - j * dis, 1.0]),
                    Vec3d([-1 + i * dis, -1, 1 - (j + 1) * dis, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, -1, 1 - (j + 1) * dis, 1.0]),
                    Vec3d([-1 + (i + 1) * dis, -1, 1 - j * dis, 1.0])
                ]


                # Extend all the squares into the array
                arrSqu.extend([Object3D(squVer1), Object3D(squVer2), Object3D(squVer3), Object3D(squVer4), Object3D(squVer5),Object3D(squVer6)])

        return arrSqu

    @staticmethod
    def catmullClark(subFaceArr,a):         #subFaceArr holds the every face and their vertices info
        from object3d import Object3D
        if a==0:
            for i in range(len(subFaceArr)):
                for j in range(len(subFaceArr[i])):
                    subFaceArr[i][j]=Vec3d(subFaceArr[i][j])
                subFaceArr[i]=Object3D(subFaceArr[i])


            return subFaceArr
        else:
            return Geo.pointMaker(subFaceArr,a-1)



    @staticmethod
    def pointMaker(subFacesArr,a):  #subFacesArr holds the original points
        subEdges = []
        subNewFacesArr=[]
        subFaces = []
        def faceCalc(faceArr):
            xVal, yVal, zVal = 0, 0, 0
            xVal += faceArr[0][0] + faceArr[1][0] + faceArr[2][0] + faceArr[3][0]
            yVal += faceArr[0][1] + faceArr[1][1] + faceArr[2][1] + faceArr[3][1]
            zVal += faceArr[0][2] + faceArr[1][2] + faceArr[2][2] + faceArr[3][2]
            return [xVal / 4, yVal / 4, zVal / 4]

        def edgeCalc(edgeArr):
            xVal, yVal, zVal = 0, 0, 0

            xVal += edgeArr[0][0] + edgeArr[1][0]
            yVal += edgeArr[0][1] + edgeArr[1][1]
            zVal += edgeArr[0][2] + edgeArr[1][2]
            return [xVal / 2, yVal / 2, zVal / 2]

        for i in range(len(subFacesArr)):
            v0 = subFacesArr[i][0]
            v1 = subFacesArr[i][1]
            v2 = subFacesArr[i][2]
            v3 = subFacesArr[i][3]


            vMid = faceCalc([v0,v1,v2,v3])
            edge0 = edgeCalc(sorted([v0, v1]))
            edge1 = edgeCalc(sorted([v1, v2]))
            edge2 = edgeCalc(sorted([v2, v3]))
            edge3 = edgeCalc(sorted([v3, v0]))



            subFaces.append([v0,v1,v2,v3,vMid])
            subEdges.extend([[edge0, [v0, v1]], [edge1, [v1, v2]],[edge2, [v2, v3]], [edge3, [v3, v0]]])

            f0 = [v0, edge0, vMid, edge3]
            f1 = [v1,edge1,vMid,edge0]
            f2 = [v2,edge2,vMid,edge1]
            f3 = [v3,edge3,vMid,edge2]
            subNewFacesArr.append(f0)
            subNewFacesArr.append(f1)
            subNewFacesArr.append(f2)
            subNewFacesArr.append(f3)

        for i in range(len(subEdges)):
            e0=subEdges[i][1][0]
            e1=subEdges[i][1][1]
            for face in subFaces:
                vertices = face[:4]
                if e0 in vertices and e1 in vertices:
                    face_point = face[4]
                    subEdges[i].append(face_point)

        for i in range(len(subNewFacesArr)):
            e1 = subNewFacesArr[i][1]
            e3 = subNewFacesArr[i][3]
            face1 = [0, 0, 0]
            face3 = [0, 0, 0]

            for j in range(len(subEdges)):
                if e1 == subEdges[j][0]:
                    face1 = [face1[k] + subEdges[j][2][k] for k in range(3)]

                if e3 == subEdges[j][0]:
                    face3 = [face3[k] + subEdges[j][2][k] for k in range(3)]

            e1_updated = [(face1[k] + 2*e1[k]) / 4 for k in range(3)]
            e3_updated = [(face3[k] + 2*e3[k]) / 4 for k in range(3)]

            subNewFacesArr[i][1] = e1_updated
            subNewFacesArr[i][3] = e3_updated

        for i in range(len(subNewFacesArr)):
            n=0
            edge=[0,0,0]
            face = [0, 0, 0]
            org_point = subNewFacesArr[i][0]
            for j in range(len(subFaces)):
                for k in range(4):
                    if subFaces[j][k]==org_point:
                        face = [face[k] + subFaces[j][4][k] for k in range(3)]
                        #print(face)
                        n+=1
                        #print(n)
            face = [coord / n for coord in face]

            for j in range(len(subEdges)):
                if subEdges[j][1][0]==org_point or subEdges[j][1][1]==org_point:
                    edge = [edge[k] + subEdges[j][0][k] for k in range(3)]

            edge = [coord / n for coord in edge]

               #this is valence now

            x = (face[0] + edge[0] + (n - 3) * org_point[0]) / n
            y = (face[1] + edge[1] + (n - 3) * org_point[1]) / n
            z = (face[2] + edge[2] + (n - 3) * org_point[2]) / n
            subNewFacesArr[i][0] = [x, y, z]


        return Geo.catmullClark(subNewFacesArr,a)


    @staticmethod
    def parser(faces,a=1):


        subFacesArr = []


        for face in faces:
            #generic quad parser
            #work logic: It calculates what should be the next point based on the distance diffrence to next vertex respect j and k (j is the vertex which is not adjacent to our vertex and k is the next one counter clockwise)


            f0 = faces[face[0]]
            f1 = faces[face[1]]
            f2 = faces[face[2]]
            f3 = faces[face[3]]

            for j in range(2**a):
                for k in range(2**a):
                    s0 = k / 2**a#**a
                    s1 = (k + 1) / 2**a#**a
                    t0 = j / 2**a#**a
                    t1 = (j + 1) / 2**a#**a

                    v0 = Vec3d.interpolate_bilinear(Vec3d(f0), Vec3d(f1), Vec3d(f2), Vec3d(f3), s0, t0)
                    v1 = Vec3d.interpolate_bilinear(Vec3d(f0), Vec3d(f1), Vec3d(f2), Vec3d(f3), s1, t0)
                    v2 = Vec3d.interpolate_bilinear(Vec3d(f0), Vec3d(f1), Vec3d(f2), Vec3d(f3), s1, t1)
                    v3 = Vec3d.interpolate_bilinear(Vec3d(f0), Vec3d(f1), Vec3d(f2), Vec3d(f3), s0, t1)

                    #print(v0.x())
                    new_face=[v0, v1, v2, v3]

                    subFacesArr.append(new_face)#subdivided objects faces and faces vertexes                   Object3D([v0, v1, v2, v3])
        return subFacesArr

    @staticmethod
    def makePyramid(a):
        from object3d import Object3D
        arrPyr=[]
        # Base face (y = 0)
        bottom = [
            Vec3d([-1.0, 0.0, 1.0, 1.0]),  # Front-left
            Vec3d([1.0, 0.0, 1.0, 1.0]),  # Front-right
            Vec3d([1.0, 0.0, -1.0, 1.0]),  # Back-right
            Vec3d([-1.0, 0.0, -1.0, 1.0])  # Back-left
        ]
        arrPyr.append(Object3D(bottom))

        dis = 1 / 2 ** a
        # up pointing triangles
        for j in range(0, 2 ** a):
            for i in range(0, j + 1):
                v1 = [
                    Vec3d([0 - (j + 1) * dis + 2 * (i + 1) * dis, 1 - dis * (j + 1), 0 + (j + 1) * dis, 1]),  # right
                    Vec3d([0 - j * dis + 2 * i * dis, 1 - dis * j, 0 + j * dis, 1]),
                    Vec3d([0 - (j + 1) * dis + 2 * i * dis, 1 - dis * (j + 1), 0 + (j + 1) * dis, 1])  # left
                ]
                v2 = [
                    Vec3d([0 - (j + 1) * dis + 2 * (i + 1) * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis, 1]),  # right
                    Vec3d([0 - j * dis + 2 * i * dis, 1 - dis * j, 0 - j * dis, 1]),
                    Vec3d([0 - (j + 1) * dis + 2 * i * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis, 1])  # left
                ]
                v3 = [
                    Vec3d([0 + (j + 1) * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis + 2 * (i + 1) * dis, 1]),  # right
                    Vec3d([0 + j * dis, 1 - dis * j, 0 - j * dis + 2 * i * dis, 1]),
                    Vec3d([0 + (j + 1) * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis + 2 * i * dis, 1])  # left
                ]
                v4 = [
                    Vec3d([0 - (j + 1) * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis + 2 * (i + 1) * dis, 1]),  # right
                    Vec3d([0 - j * dis, 1 - dis * j, 0 - j * dis + 2 * i * dis, 1]),
                    Vec3d([0 - (j + 1) * dis, 1 - dis * (j + 1), 0 - (j + 1) * dis + 2 * i * dis, 1])  # left
                ]
                arrPyr.extend([Object3D(v1), Object3D(v2), Object3D(v3), Object3D(v4)])
        # down pointing triangles
        for j in range(0, 2 ** a):
            for i in range(0, j):
                v1 = [
                    Vec3d([0 + 2 * i * dis - j * dis, 1 - j * dis, 0 + j * dis, 1]),  # left
                    Vec3d([0 + 2 * (i + 1) * dis - j * dis, 1 - j * dis, 0 + j * dis, 1]),  # right
                    Vec3d([0 + 2 * i * dis - (j - 1) * dis, 1 - (j + 1) * dis, 0 + (j + 1) * dis, 1])  # bottom
                ]
                v2 = [
                    Vec3d([0 + 2 * i * dis - j * dis, 1 - j * dis, 0 - j * dis, 1]),  # left
                    Vec3d([0 + 2 * (i + 1) * dis - j * dis, 1 - j * dis, 0 - j * dis, 1]),  # right
                    Vec3d([0 + 2 * i * dis - (j - 1) * dis, 1 - (j + 1) * dis, 0 - (j + 1) * dis, 1])  # bottom
                ]
                v3 = [
                    Vec3d([0 + j * dis, 1 - j * dis, 0 + 2 * i * dis - j * dis, 1]),  # left
                    Vec3d([0 + j * dis, 1 - j * dis, 0 + 2 * (i + 1) * dis - j * dis, 1]),  # right
                    Vec3d([0 + (j + 1) * dis, 1 - (j + 1) * dis, 0 + 2 * i * dis - (j - 1) * dis, 1])  # bottom
                ]
                v4 = [
                    Vec3d([0 - j * dis, 1 - j * dis, 0 + 2 * i * dis - j * dis, 1]),  # left
                    Vec3d([0 - j * dis, 1 - j * dis, 0 + 2 * (i + 1) * dis - j * dis, 1]),  # right
                    Vec3d([0 - (j + 1) * dis, 1 - (j + 1) * dis, 0 + 2 * i * dis - (j - 1) * dis, 1])  # bottom
                ]
                arrPyr.extend([Object3D(v1), Object3D(v2), Object3D(v3), Object3D(v4)])
        return arrPyr



    @staticmethod
    def makePlane(a):
        from object3d import Object3D
        arrPlane = []

        dis = (2 / (2 ** a))

        for j in range(0, 2 ** a):
            for i in range(0, 2 ** a):
                planeVer = [
                    Vec3d([-1 + i * dis, 1 - j * dis, 0, 1.0]),  # top left
                    Vec3d([-1 + i * dis, 1 - (j + 1) * dis, 0, 1.0]),  # bot left
                    Vec3d([-1 + (i + 1) * dis, 1 - (j + 1) * dis, 0, 1.0]),  # bot rigth
                    Vec3d([-1 + (i + 1) * dis, 1 - j * dis, 0, 1.0])  # top right
                ]
                v1 = Object3D(planeVer)
                arrPlane.append(v1)

        return arrPlane



    @staticmethod
    def makeSphere(a):
        from object3d import Object3D
        arrSphere = []

        if a<0:
            a=0

        stacks = 2**a
        slices = 2**a
        d_theta = pi / stacks
        d_phi = pi / slices
        #x=sin(m)*cos(n)
        #y=sin(m)*sin(n)
        #z=cos(m)
        for i in range(stacks):
            theta1 = i * d_theta
            theta2 = (i + 1) * d_theta

            for j in range(slices):
                phi1 = j * d_phi
                phi2 = (j + 1) * d_phi

                # Vertices for the first hemisphere
                v1 = Vec3d([sin(theta1) * cos(phi1), sin(theta1) * sin(phi1), cos(theta1), 1])
                v2 = Vec3d([sin(theta2) * cos(phi1), sin(theta2) * sin(phi1), cos(theta2), 1])
                v3 = Vec3d([sin(theta2) * cos(phi2), sin(theta2) * sin(phi2), cos(theta2), 1])
                v4 = Vec3d([sin(theta1) * cos(phi2), sin(theta1) * sin(phi2), cos(theta1), 1])

                # Add the triangles (v1, v2, v3) and (v1, v3, v4)
                arrSphere.append(Object3D([v1, v2, v3]))
                arrSphere.append(Object3D([v1, v3, v4]))

        for i in range(stacks):
            theta1 = i * d_theta + pi
            theta2 = (i + 1) * d_theta + pi

            for j in range(slices):
                phi1 = j * d_phi
                phi2 = (j + 1) * d_phi


                v1 = Vec3d([sin(theta1) * cos(phi1), sin(theta1) * sin(phi1), cos(theta1), 1])
                v2 = Vec3d([sin(theta2) * cos(phi1), sin(theta2) * sin(phi1), cos(theta2), 1])
                v3 = Vec3d([sin(theta2) * cos(phi2), sin(theta2) * sin(phi2), cos(theta2), 1])
                v4 = Vec3d([sin(theta1) * cos(phi2), sin(theta1) * sin(phi2), cos(theta1), 1])


                arrSphere.append(Object3D([v1, v2, v3]))
                arrSphere.append(Object3D([v1, v3, v4]))

        return arrSphere

    @staticmethod
    def makeSphere2(a):
        from object3d import Object3D
        arrSphere = []

        #x=sin(m)*cos(n)
        #y=sin(m)*sin(n)
        #z=cos(m)
        dn=pi/2/2**a+7   #this is the angle we will use      -   n is the x y angle
        n=pi/2
        for j in range(1, 2 ** a+1):
            m=pi/2
            for i in range(0, j + 1):#d1 is for the upper point of the triangle
                if j==0:
                    dm1=0
                else:
                    dm1=m/j
                dm2=m/(j+1)      #d2 is for the bottom of the triangle
                v1 = [
                    Vec3d([sin(m-dm1*i)*cos(n), sin(m-dm1*i)*sin(n), cos(m-dm1*i), 1]),
                    Vec3d([sin(m-dm2*i)*cos(n-dn), sin(m-dm2*i)*sin(n-dn), cos(m-dm2*i), 1]),
                    Vec3d([sin(m-dm2*(i+1))*cos(n-dn), sin(m-dm2*(i+1))*sin(n-dn), cos(m-dm2*(i+1)), 1])
                ]
                arrSphere.extend([Object3D(v1)])
            n=n-dn


        dn=pi/2/2**a+7   #this is the angle we will use      -   n is the x y angle
        n=pi
        for j in range(1, 2 ** a+1):
            m=pi/2
            for i in range(0, j + 1):#d1 is for the upper point of the triangle
                if j==0:
                    dm1=0
                else:
                    dm1=m/j
                dm2=m/(j+1)      #d2 is for the bottom of the triangle
                v1 = [#this is a face
                    Vec3d([sin(m-dm1*i)*cos(n), sin(m-dm1*i)*sin(n), -cos(m-dm1*i), 1]),
                    Vec3d([sin(m-dm2*i)*cos(n-dn), sin(m-dm2*i)*sin(n-dn), -cos(m-dm2*i), 1]),
                    Vec3d([sin(m-dm2*(i+1))*cos(n-dn), sin(m-dm2*(i+1))*sin(n-dn), -cos(m-dm2*(i+1)), 1])
                ]
                arrSphere.extend([Object3D(v1)])
            n=n-dn

        return arrSphere

    @staticmethod
    def drawObject(arrObj,lines):
        for i in range(len(arrObj)):
            face = arrObj[i]
            vertices = face.getVertices()
            num_vertices = len(vertices)


            if num_vertices == 3:
                glBegin(GL_TRIANGLES)
            elif num_vertices == 4:
                glBegin(GL_QUADS)
            else:
                glBegin(GL_POLYGON)

            r = ((i + 2.5) % 3) / 3.0
            g = ((i + 6.5) % 7) / 7.0
            b = ((i + 4.5) % 5) / 5.0

            glColor3f(r, g, b)

            for vertex in vertices:
                transformed = cam.returnViewMatrix().matrixVecMult(vertex)
                glVertex3f(transformed.x(), transformed.y(), transformed.z())
            glEnd()

            if lines:
                glColor3f(1.0, 1.0, 1.0)
                glBegin(GL_LINE_LOOP if num_vertices == 3 else GL_LINE_STRIP)
                for vertex in vertices:
                    transformed = cam.returnViewMatrix().matrixVecMult(vertex)
                    glVertex3f(transformed.x(), transformed.y(), transformed.z())

                if num_vertices == 4:
                    transformed = cam.returnViewMatrix().matrixVecMult(vertices[0])
                    glVertex3f(transformed.x(), transformed.y(), transformed.z())
                glEnd()

    @staticmethod
    def resetGeoCam(x,y,z):
        cam.changeCamPos(x, y, z)
        cam.lookInto(0.0,0.0,-8.0)
    @staticmethod
    def rotateAround(x,y):
        cam.rotateAround(x,y)
    @staticmethod
    def dragZ(z):
        cam.dragZ(z)
    @staticmethod
    def changeLookInto(x,y,z):
        cam.changeLookInto(x,y,z)

