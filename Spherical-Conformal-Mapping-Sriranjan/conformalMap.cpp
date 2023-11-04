#include <fstream>
#include <numeric>
#include <vector>
#include <math.h>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "Point.h"
#include "CTrait.h"

using namespace MeshLib;
using namespace std;
using namespace MeshLib;

#define DELTA_T 1e-2
#define DELTA_E 1e-5
#define HARMONIC 1
#define TUETTE 0

class Harmonic {

public:
	Harmonic(Solid *conMesh) {
		normMesh = conMesh;
	};

	void Star();
	void harmonicMapConjugate();
	void harmonicMapGD();
	void Tuette();
	void set_up();

	std::vector<double> computeVertexAngle();

protected:
	double energyCompute(int type);
	void gradientCompute(int type);
	void gradientConjugate();
	void setNormal();
	void kuv();
	void meshUpdate();
	void updateMassCenter();

	Solid * normMesh;
};

void Harmonic::Star() {
	Point center(0, 0, 0);
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		center += vertex->point();
	}
	center /= normMesh->numVertices();
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
		v_n(vertex) = vertex->point();
	}
}

double Harmonic::energyCompute(int type) {
	double energy = 0;
	for (SolidEdgeIterator eiter(normMesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Solid::tVertex v1 = normMesh->edgeVertex1(edge);
		Solid::tVertex v2 = normMesh->edgeVertex2(edge);
		Point uv = v1->point() - v2->point();
		if (type == HARMONIC) {
			energy += e_k(edge) * uv.norm2();
		}
		else energy += uv.norm2();
	}
	return energy;
}

void Harmonic::gradientCompute(int type) {
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex sv = *viter;
		Point gradient = Point(0, 0, 0);
		for (VertexVertexIterator vviter(sv); !vviter.end(); ++vviter) {
			Vertex * tv = *vviter;
			if (type == HARMONIC)
				gradient += (sv->point() - tv->point()) * e_k(normMesh->vertexEdge(sv, tv));
			else gradient += (sv->point() - tv->point());
		}
		v_dv(sv) = gradient;
		v_abdv(sv) = v_dv(sv) - v_n(sv) * (v_dv(sv) * v_n(sv));
	}
}

void Harmonic::gradientConjugate() {
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		v_beta(vertex) = v_abdv(vertex)*v_abdv(vertex) / (v_fabdv(vertex)*v_fabdv(vertex));
		v_s(vertex) = v_s(vertex)*v_beta(vertex) - v_abdv(vertex);
	}
	double alpha = 1e-6;
	double e0 = energyCompute(HARMONIC);
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v_mp(v) = v->point();
	}
	while (alpha<1e-2) {
		for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
			Solid::tVertex v = *viter;
			v->point() = v->point() + v_s(v)*alpha;
			v->point() /= v->point().norm();
		}
		double nenergy = energyCompute(HARMONIC);
		if (nenergy < e0) {
			for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
				Solid::tVertex v = *viter;
				v_mp(v) = v->point();
			}
		}
		alpha *= 2;
	}
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v_mp(v);
	}
}

void Harmonic::updateMassCenter() {
	for (SolidFaceIterator fiter(normMesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;

		Point p[3];
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex * v = *fviter;
			p[i++] = v->point();
		}

		Point n = (p[1] - p[0]) ^ (p[2] - p[0]);
		f_a(face) = n.norm() / 2.0;
	}
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		double area = 0;
		for (VertexFaceIterator vfiter(vertex); !vfiter.end(); ++vfiter) {
			Face* vf = *vfiter;
			area += f_a(vf);
		}
		v_a(vertex) = area / 3.0;
	}

	Point center(0, 0, 0);
	double mass = 0;
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		center += vertex->point() * v_a(vertex);
		mass += v_a(vertex);
	}
	center /= mass;

	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
	}
}

void Harmonic::meshUpdate() {
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v->point()- v_abdv(v)*DELTA_T;
		v->point() /= v->point().norm();
	}
}

void Harmonic::harmonicMapConjugate() {
	double te = energyCompute(HARMONIC);
	double oe = 1000;
	gradientCompute(HARMONIC);
	meshUpdate();
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		v_fabdv(vertex) = v_abdv(vertex);
		v_s(vertex) = v_abdv(vertex);
	}
	double he = energyCompute(HARMONIC);
	int i = 0;
	while (oe - he > DELTA_E) {
		gradientCompute(HARMONIC);
		gradientConjugate();
		oe = he;
		he = energyCompute(HARMONIC);
		if (i++ < 10 or i % 100 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::harmonicMapGD() {
	double he = energyCompute(HARMONIC);
	double oe = 1000;
	int i = 0;
	while (oe - he > DELTA_E) {
		gradientCompute(HARMONIC);
		meshUpdate();
		updateMassCenter();
		oe = he;
		he = energyCompute(HARMONIC);
		if (i++ < 10 or (i < 500 and i % 100 == 0) or i % 1000 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::Tuette() {
	double te = energyCompute(TUETTE);
	double oe = 1000;
	int i = 0;
	while (oe-te > DELTA_E) {
		gradientCompute(TUETTE);
		meshUpdate();
		oe = te;
		te = energyCompute(TUETTE);
		if (i++ < 10 or (i < 500 and i % 100 == 0) or i % 1000 == 0) {
			std::cout << i << ", energy: " << te << ", diff: "<< oe - te << std::endl;	
		}
	}
}

void Harmonic::setNormal() {
	for (SolidVertexIterator viter(normMesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->trait() = new CVertexTrait;
	}
	for (SolidFaceIterator fiter(normMesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;
		face->trait() = new CFaceTrait;
		Point fp[3]; 
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex* fv = *fviter;
			fp[i++] = fv->point();
		}
		Point normal = (fp[1] - fp[0]) ^ (fp[2] - fp[0]);
		f_n(face) = normal / normal.norm();
	}
}

void Harmonic::kuv() {
	for (SolidEdgeIterator eiter(normMesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		edge->trait() = new CEdgeTrait;
		e_k(edge) = 1.0;
	}

	double sum = 0;
	for (SolidEdgeIterator eiter(normMesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Point p1 = normMesh->edgeVertex1(edge)->point();
		Point p2 = normMesh->edgeVertex2(edge)->point();
		Point p3 = edge->halfedge(0)->he_next()->target()->point();

		double alpha, beta;
		alpha = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		p3 = edge->halfedge(0)->he_sym()->he_next()->target()->point();
		beta = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		e_k(edge) = alpha + beta;
		sum += e_k(edge);
	}
}

void Harmonic::set_up() {
	setNormal();
	kuv();
}

std::vector<double> Harmonic::computeVertexAngle() {
	std::vector<double> angles;
	for (SolidFaceIterator fiter(normMesh); !fiter.end(); ++fiter) {
		Solid::tFace f = *fiter;
		Point p[3];
		int i = 0;
		for (FaceVertexIterator fviter(f); !fviter.end(); ++fviter) {
			Vertex * v = *fviter;
			p[i++] = v->point();
		}
		Point p1 = p[0] - p[1];
		Point p2 = p[1] - p[2];
		Point p3 = p[2] - p[0];

		angles.push_back(acos(p1 * (-p3) / (p1.norm() * (-p3).norm())));
		angles.push_back(acos((-p1) * p2 / ((-p1).norm() * p2.norm())));
		angles.push_back(acos(p3 * (-p2) / (p3.norm() * (-p2).norm())));
	}
	return angles;
}

int main(int argc, char *argv[]) {
	Solid nmesh;
	OBJFileReader of;
	ifstream in(argv[1]);
	of.readToSolid(&nmesh, in);

	Harmonic shm(&nmesh);
	cout << "--> Setting up" << endl;
	vector<double> angles_before = shm.computeVertexAngle();
	shm.set_up();
	
	shm.Star();
	cout << "--> Computing Tuette map..." << endl;
	
	shm.Tuette();
	cout << "--> Computing Harmonic map..." << endl;
	
	shm.harmonicMapConjugate();
	string out = string(argv[2]) + "_harmonic.obj";
	
	of.writeToObj(&nmesh, out);

	return 0;
}