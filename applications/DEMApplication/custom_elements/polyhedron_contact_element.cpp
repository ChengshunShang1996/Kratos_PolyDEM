/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

// Project includes
#include "custom_elements/polyhedron_contact_element.h"
#include "utilities/math_utils.h"
#include "DEM_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PolyhedronContactElement::PolyhedronContactElement(IndexType NewId, 
                                                    PolyhedronParticle* PolyhedronParticle1, 
                                                    PolyhedronParticle* PolyhedronParticle2): 
mId(NewId),
mPolyhedronParticle1(PolyhedronParticle1),
mPolyhedronParticle2(PolyhedronParticle2)
{

}

PolyhedronContactElement::PolyhedronContactElement(IndexType NewId, 
                                                    PolyhedronParticle* PolyhedronParticle1, 
                                                    PolyhedronParticle* PolyhedronParticle2, 
                                                    PropertiesType::Pointer pProperties):
mId(NewId),
mPolyhedronParticle1(PolyhedronParticle1),
mPolyhedronParticle2(PolyhedronParticle2),
mpProperties(pProperties)
{

}

//create contact elements instances.
PolyhedronContactElement::Pointer PolyhedronContactElement::Create(IndexType NewId, PolyhedronParticle* PolyhedronParticle1, PolyhedronParticle* PolyhedronParticle2, PropertiesType::Pointer pProperties) const
{
    return PolyhedronContactElement::Pointer ( new PolyhedronContactElement(NewId, PolyhedronParticle1, PolyhedronParticle2, pProperties));
}

PolyhedronContactElement::~PolyhedronContactElement()
{
}

void PolyhedronContactElement::Initialize(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    mFailureState = 0.0;
    mContactForce[0] = 0.0;
    mContactForce[1] = 0.0;
    mContactForce[2] = 0.0;
    mRotationalMoment[0] = 0.0;
    mRotationalMoment[1] = 0.0;
    mRotationalMoment[2] = 0.0;
	SetDeleteFlag(false);
	mGlobalDamping = r_process_info[GLOBAL_DAMPING];
	mTangentialElasticContactForce = Vector3(0.0, 0.0, 0.0);

    /*
    array_1d<double, 3> vector_of_zeros(3, 0.0);
    this->SetValue(LOCAL_CONTACT_FORCE, vector_of_zeros);
    this->SetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT, vector_of_zeros);
    this->SetValue(FAILURE_CRITERION_STATE, 0.0);
    this->SetValue(CONTACT_RADIUS, 0.0);*/

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::PrepareForPrinting() {
    KRATOS_TRY

    /*
    this->GetValue(LOCAL_CONTACT_FORCE)[0] = mContactForce[0];
    this->GetValue(LOCAL_CONTACT_FORCE)[1] = mContactForce[1];
    this->GetValue(LOCAL_CONTACT_FORCE)[2] = mContactForce[2];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[0] = mContactForce[0];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[1] = mRotationalMoment[1];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[2] = mRotationalMoment[2];
    this->GetValue(FAILURE_CRITERION_STATE)= mFailureCriterionState;
    this->GetValue(CONTACT_RADIUS)         = mContactRadius;*/

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::InitializeSolutionStep(const ProcessInfo& r_process_info )
{
    //
}

void PolyhedronContactElement::FinalizeSolutionStep(const ProcessInfo& r_process_info) {

}

void PolyhedronContactElement::CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) {
    KRATOS_TRY

    if (GJK()) {
        
		auto& central_node_1 = mPolyhedronParticle1->GetGeometry()[0];
		auto& central_node_2 = mPolyhedronParticle2->GetGeometry()[0];
		Vector3 coll1Pos = {central_node_1[0], central_node_1[1], central_node_1[2]};
		Vector3 coll2Pos = {central_node_2[0], central_node_2[1], central_node_2[2]};
		Vector3 contact_m = (mContactPoint1 + mContactPoint2)/2;

		Vector3 contact_force(0.0, 0.0, 0.0);
		Vector3 contact_moment_1(0.0, 0.0, 0.0);
		Vector3 contact_moment_2(0.0, 0.0, 0.0);

		//double kn = 100000.0;
		//contact_force = mOverlapVector * kn;

		ClonePolyhedronDiscontinuumConstitutiveLawWithNeighbour();
    	mPolyhedronDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, mPolyhedronParticle1, mPolyhedronParticle2, 
																mOverlapVector, contact_force, mTangentialElasticContactForce);

		Vector3 torque_arm_1 = contact_m - coll1Pos;
		Vector3 torque_arm_2 = contact_m - coll2Pos;
		contact_moment_1 = Vector3::Cross(torque_arm_1, contact_force);
		contact_moment_2 = Vector3::Cross(torque_arm_2, -contact_force);

		array_1d<double,3>& total_forces_1 = central_node_1.FastGetSolutionStepValue(TOTAL_FORCES);
		array_1d<double,3>& total_moment_1 = central_node_1.FastGetSolutionStepValue(PARTICLE_MOMENT);

		total_forces_1[0] = contact_force[0];
		total_forces_1[1] = contact_force[1];
		total_forces_1[2] = contact_force[2];

		total_moment_1[0] = contact_moment_1[0];
		total_moment_1[1] = contact_moment_1[1];
		total_moment_1[2] = contact_moment_1[2];

		array_1d<double,3>& total_forces_2 = central_node_2.FastGetSolutionStepValue(TOTAL_FORCES);
		array_1d<double,3>& total_moment_2 = central_node_2.FastGetSolutionStepValue(PARTICLE_MOMENT);

		total_forces_2[0] = -contact_force[0];
		total_forces_2[1] = -contact_force[1];
		total_forces_2[2] = -contact_force[2];

		total_moment_2[0] = contact_moment_2[0];
		total_moment_2[1] = contact_moment_2[1];
		total_moment_2[2] = contact_moment_2[2];

		/*
		total_forces[0] = contact_force[0] + additional_forces[0];
		total_forces[1] = contact_force[1] + additional_forces[1];
		total_forces[2] = contact_force[2] + additional_forces[2];

		total_moment[0] = mContactMoment[0] + additionally_applied_moment[0];
		total_moment[1] = mContactMoment[1] + additionally_applied_moment[1];
		total_moment[2] = mContactMoment[2] + additionally_applied_moment[2];*/

		ApplyGlobalDampingToContactForcesAndMoments(mPolyhedronParticle1, total_forces_1, total_moment_1);
		ApplyGlobalDampingToContactForcesAndMoments(mPolyhedronParticle2, total_forces_2, total_moment_2);

		#ifdef KRATOS_DEBUG
		DemDebugFunctions::CheckIfNan(total_forces_1, "NAN in Total Forces in RHS of Particle");
		DemDebugFunctions::CheckIfNan(total_forces_2, "NAN in Total Torque in RHS of Particle");
		DemDebugFunctions::CheckIfNan(total_moment_1, "NAN in Total Forces in RHS of Particle");
		DemDebugFunctions::CheckIfNan(total_moment_2, "NAN in Total Torque in RHS of Particle");
		#endif
    }

    KRATOS_CATCH( "" )
}

bool PolyhedronContactElement::GJK()
{
	KRATOS_TRY

	auto& central_node_1 = mPolyhedronParticle1->GetGeometry()[0];
    auto& central_node_2 = mPolyhedronParticle2->GetGeometry()[0];

    Vector3 coll1Pos = {central_node_1[0], central_node_1[1], central_node_1[2]};
    Vector3 coll2Pos = {central_node_2[0], central_node_2[1], central_node_2[2]};

	Vector3 search_dir = coll1Pos - coll2Pos; //initial search direction between colliders
	Point a, b, c, d; //Simplex: just a set of points (a is always most recently added)

	//Get initial point for simplex
	//Point c;
	CalculateSearchPoint(c, search_dir);
	search_dir = -c.p; //search in direction of origin

	//Get second point for a line segment simplex
	//Point b;
	CalculateSearchPoint(b, search_dir);

	if (Vector3::Dot(b.p, search_dir) < 0) {
		return false;
	}//we didn't reach the origin, won't enclose it

	search_dir = Vector3::Cross(Vector3::Cross(c.p - b.p, -b.p), c.p - b.p); //search perpendicular to line segment towards origin
	if (search_dir == Vector3(0.0, 0.0, 0.0)) { //origin is on this line segment
		//Apparently any normal search vector will do?
		search_dir = Vector3::Cross(c.p - b.p, Vector3(1.0, 0.0, 0.0)); //normal with x-axis
		if (search_dir == Vector3(0.0, 0.0, 0.0)){
			search_dir = Vector3::Cross(c.p - b.p, Vector3(0.0, 0.0, -1.0)); //normal with z-axis
		}
	}
			
	int simp_dim = 2; //simplex dimension

	for (int iterations = 0; iterations < GJK_MAX_NUM_ITERATIONS; iterations++)
	{
		//Point a;
		CalculateSearchPoint(a, search_dir);

		if (Vector3::Dot(a.p, search_dir) < 0) {
			return false;
		}//we didn't reach the origin, won't enclose it

		simp_dim++;
		if (simp_dim == 3) {
			update_simplex3(a, b, c, d, simp_dim, search_dir);
		}
		else if (update_simplex4(a, b, c, d, simp_dim, search_dir)) {
			EPA(a, b, c, d);
			return true;
		}
	}//endfor

	return false;

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::update_simplex3(Point& a, Point& b, Point& c, Point& d, int& simp_dim, Vector3& search_dir)
{
	KRATOS_TRY

    Vector3 n = Vector3::Cross(b.p - a.p, c.p - a.p); //triangle's normal
	Vector3 AO = -a.p; //direction to origin

	//Determine which feature is closest to origin, make that the new simplex
	simp_dim = 2;
	if (Vector3::Dot(Vector3::Cross(b.p - a.p, n), AO) > 0) { //Closest to edge AB
		c = a;
		//simp_dim = 2;
		search_dir = Vector3::Cross(Vector3::Cross(b.p - a.p, AO), b.p - a.p);
		return;
	}
	if (Vector3::Dot(Vector3::Cross(n, c.p - a.p), AO) > 0) { //Closest to edge AC
		b = a;
		//simp_dim = 2;
		search_dir = Vector3::Cross(Vector3::Cross(c.p - a.p, AO), c.p - a.p);
		return;
	}

	simp_dim = 3;
	if (Vector3::Dot(n, AO) > 0) { //Above triangle
		d = c;
		c = b;
		b = a;
		search_dir = n;
		return;
	}
	//else //Below triangle
	d = b;
	b = a;
	search_dir = -n;
	return;

    KRATOS_CATCH( "" )
}

bool PolyhedronContactElement::update_simplex4(Point& a, Point& b, Point& c, Point& d, int& simp_dim, Vector3& search_dir)
{
	KRATOS_TRY
    // a is peak/tip of pyramid, BCD is the base (counterclockwise winding order)
	//We know a priori that origin is above BCD and below a

	//Get normals of three new faces
	Vector3 ABC = Vector3::Cross(b.p - a.p, c.p - a.p);
	Vector3 ACD = Vector3::Cross(c.p - a.p, d.p - a.p);
	Vector3 ADB = Vector3::Cross(d.p - a.p, b.p - a.p);

	Vector3 AO = -a.p; //dir to origin
	simp_dim = 3; //hoisting this just cause

	if (Vector3::Dot(ABC, AO) > 0) { //In front of ABC
		d = c;
		c = b;
		b = a;
		search_dir = ABC;
		return false;
	}

	if (Vector3::Dot(ACD, AO) > 0) { //In front of ACD
		b = a;
		search_dir = ACD;
		return false;
	}
	if (Vector3::Dot(ADB, AO) > 0) { //In front of ADB
		c = d;
		d = b;
		b = a;
		search_dir = ADB;
		return false;
	}

	//else inside tetrahedron; enclosed!
	return true;

    KRATOS_CATCH( "" )
}

//Expanding Polytope Algorithm
void PolyhedronContactElement::EPA(Point& a, Point& b, Point& c, Point& d)
{
	KRATOS_TRY

    Point faces[EPA_MAX_NUM_FACES][4]; //Array of faces, each with 3 verts and a normal

	Vector3 VertexA[3];
	Vector3 VertexB[3];

	//Init with final simplex from GJK
	faces[0][0] = a;
	faces[0][1] = b;
	faces[0][2] = c;
	faces[0][3].p = (Vector3::Cross(b.p - a.p, c.p - a.p)).Normalised(); //ABC
	faces[1][0] = a;
	faces[1][1] = c;
	faces[1][2] = d;
	faces[1][3].p = (Vector3::Cross(c.p - a.p, d.p - a.p)).Normalised(); //ACD
	faces[2][0] = a;
	faces[2][1] = d;
	faces[2][2] = b;
	faces[2][3].p = (Vector3::Cross(d.p - a.p, b.p - a.p)).Normalised(); //ADB
	faces[3][0] = b;
	faces[3][1] = d;
	faces[3][2] = c;
	faces[3][3].p = (Vector3::Cross(d.p - b.p, c.p - b.p)).Normalised(); //BDC

	int num_faces = 4;
	int closest_face;

	for (int iterations = 0; iterations < EPA_MAX_NUM_ITERATIONS; iterations++) {
		//Find face that's closest to origin
		double min_dist = Vector3::Dot(faces[0][0].p, faces[0][3].p);
		closest_face = 0;
		for (int i = 1; i < num_faces; i++) {
			double dist = Vector3::Dot(faces[i][0].p, faces[i][3].p);
			if (dist < min_dist) {
				min_dist = dist;
				closest_face = i;
			}
		}

		//search normal to face that's closest to origin
		Vector3 search_dir = faces[closest_face][3].p;

		Point p;
		CalculateSearchPoint(p, search_dir);

		if (Vector3::Dot(p.p, search_dir) - min_dist < EPA_TOLERANCE) {

			PolyPlane closestPlane = PolyPlane::PlaneFromTri(faces[closest_face][0].p, faces[closest_face][1].p, faces[closest_face][2].p); //plane of closest triangle face
			Vector3 projectionPoint = closestPlane.ProjectPointOntoPlane(Vector3(0, 0, 0)); //projecting the origin onto the triangle(both are in Minkowski space)
			double u, v, w;
			Barycentric(faces[closest_face][0].p, faces[closest_face][1].p, faces[closest_face][2].p, projectionPoint, u, v, w); //finding the barycentric coordinate of this projection point to the triangle

			//The contact points just have the same barycentric coordinate in their own triangles which  are composed by result coordinates of support function 
			mContactPoint1 = faces[closest_face][0].a * u + faces[closest_face][1].a * v + faces[closest_face][2].a * w;
			mContactPoint2 = faces[closest_face][0].b * u + faces[closest_face][1].b * v + faces[closest_face][2].b * w;
            mOverlapVector = mContactPoint2 - mContactPoint1;
			KRATOS_WATCH(mOverlapVector);
			//Vector3 normal = (mContactPoint1 - mContactPoint2).Normalised();
            //Vector3 contact_normal = faces[closest_face][3]
            //Vector3 contact_point_minkowski = contact_normal * Vector3::Dot(p.p, search_dir)
			return;
		}

		Point loose_edges[EPA_MAX_NUM_LOOSE_EDGES][2]; //keep track of edges we need to fix after removing faces
		int num_loose_edges = 0;

		//Find all triangles that are facing p
		for (int i = 0; i < num_faces; i++)
		{
			if (Vector3::Dot(faces[i][3].p, p.p - faces[i][0].p) > 0) //triangle i faces p, remove it
			{
				//Add removed triangle's edges to loose edge list.
				//If it's already there, remove it (both triangles it belonged to are gone)
				for (int j = 0; j < 3; j++) //Three edges per face
				{
					Point current_edge[2] = { faces[i][j], faces[i][(j + 1) % 3] };
					bool found_edge = false;
					for (int k = 0; k < num_loose_edges; k++) //Check if current edge is already in list
					{
						if (loose_edges[k][1].p == current_edge[0].p && loose_edges[k][0].p == current_edge[1].p) {
							loose_edges[k][0] = loose_edges[num_loose_edges - 1][0]; //Overwrite current edge
							loose_edges[k][1] = loose_edges[num_loose_edges - 1][1]; //with last edge in list
							num_loose_edges--;
							found_edge = true;
							k = num_loose_edges; //exit loop because edge can only be shared once
						}
					}//endfor loose_edges

					if (!found_edge) { //add current edge to list
						// assert(num_loose_edges<EPA_MAX_NUM_LOOSE_EDGES);
						if (num_loose_edges >= EPA_MAX_NUM_LOOSE_EDGES) break;
						loose_edges[num_loose_edges][0] = current_edge[0];
						loose_edges[num_loose_edges][1] = current_edge[1];
						num_loose_edges++;
					}
				}

				//Remove triangle i from list
				faces[i][0] = faces[num_faces - 1][0];
				faces[i][1] = faces[num_faces - 1][1];
				faces[i][2] = faces[num_faces - 1][2];
				faces[i][3] = faces[num_faces - 1][3];
				num_faces--;
				i--;
			}//endif p can see triangle i
		}//endfor num_faces

		//Reconstruct polytope with p added
		for (int i = 0; i < num_loose_edges; i++)
		{
			// assert(num_faces<EPA_MAX_NUM_FACES);
			if (num_faces >= EPA_MAX_NUM_FACES) break;
			faces[num_faces][0] = loose_edges[i][0];
			faces[num_faces][1] = loose_edges[i][1];
			faces[num_faces][2] = p;
			faces[num_faces][3].p = (Vector3::Cross(loose_edges[i][0].p - loose_edges[i][1].p, loose_edges[i][0].p - p.p)).Normalised();

			//Check for wrong normal to maintain CCW winding
			double bias = 0.000001; //in case dot result is only slightly < 0 (because origin is on face)
			if (Vector3::Dot(faces[num_faces][0].p, faces[num_faces][3].p) + bias < 0) {
				Point temp = faces[num_faces][0];
				faces[num_faces][0] = faces[num_faces][1];
				faces[num_faces][1] = temp;
				faces[num_faces][3].p = -faces[num_faces][3].p;
			}
			num_faces++;
		}
	} //End for iterations
	std::cout<< "EPA did not converge" << std::endl;
	//Return most recent closest point
	Vector3 search_dir = faces[closest_face][3].p;

	Point p;
	CalculateSearchPoint(p, search_dir);

	PolyPlane closestPlane = PolyPlane::PlaneFromTri(faces[closest_face][0].p, faces[closest_face][1].p, faces[closest_face][2].p);
	Vector3 projectionPoint = closestPlane.ProjectPointOntoPlane(Vector3(0, 0, 0));
	double u, v, w;
	Barycentric(faces[closest_face][0].p, faces[closest_face][1].p, faces[closest_face][2].p, projectionPoint, u, v, w);
	mContactPoint1 = faces[closest_face][0].a * u + faces[closest_face][1].a * v + faces[closest_face][2].a * w;
	mContactPoint2 = faces[closest_face][0].b * u + faces[closest_face][1].b * v + faces[closest_face][2].b * w;
    mOverlapVector = mContactPoint2 - mContactPoint1;
    /*
    Vector3 localA = faces[closest_face][0].a * u + faces[closest_face][1].a * v + faces[closest_face][2].a * w;
	Vector3 localB = faces[closest_face][0].b * u + faces[closest_face][1].b * v + faces[closest_face][2].b * w;
	double penetration = (localA - localB).Length();
	Vector3 normal = (localA - localB).Normalised();*/

	return;

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::CalculateSearchPoint(Point& point, Vector3& search_dir)
{
	KRATOS_TRY

    point.b = mPolyhedronParticle2->GetFurthestPoint(search_dir);
	point.a = mPolyhedronParticle1->GetFurthestPoint(-search_dir);
	point.p = point.b - point.a;

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::Barycentric(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& p, double& u, double& v, double& w)
{
    KRATOS_TRY
    
    Vector3 v0 = b - a, v1 = c - a, v2 = p - a;
    double d00 = Vector3::Dot(v0, v0);
    double d01 = Vector3::Dot(v0, v1);
    double d11 = Vector3::Dot(v1, v1);
    double d20 = Vector3::Dot(v2, v0);
    double d21 = Vector3::Dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::SetId(IndexType NewId) { mId = NewId;}

void PolyhedronContactElement::SetPolyElement1(PolyhedronParticle* custom_poly_element){
	KRATOS_TRY

    mPolyhedronParticle1 = custom_poly_element;

    KRATOS_CATCH( "" )
}


void PolyhedronContactElement::SetPolyElement2(PolyhedronParticle* custom_poly_element){
	KRATOS_TRY

    mPolyhedronParticle2 = custom_poly_element;

    KRATOS_CATCH( "" )
};

PolyhedronParticle* PolyhedronContactElement::GetPolyElement1(){
	KRATOS_TRY

    return mPolyhedronParticle1;

    KRATOS_CATCH( "" )
}


PolyhedronParticle* PolyhedronContactElement::GetPolyElement2(){
	KRATOS_TRY

    return mPolyhedronParticle2;

    KRATOS_CATCH( "" )
};

void PolyhedronContactElement::SetDeleteFlag(bool this_flag) { mDeleteFlag = this_flag; }
bool PolyhedronContactElement::GetDeleteFlag()               { return mDeleteFlag;}

void PolyhedronContactElement::ApplyGlobalDampingToContactForcesAndMoments(PolyhedronParticle* ThisPolyhedronParticle, array_1d<double,3>& total_forces, array_1d<double,3>& total_moment) {

        KRATOS_TRY

        auto& central_node = ThisPolyhedronParticle->GetGeometry()[0];

		const array_1d<double, 3> velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3> angular_velocity = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);

        if (central_node.IsNot(DEMFlags::FIXED_VEL_X)) {
            total_forces[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[0] * velocity[0]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Y)) {
            total_forces[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[1] * velocity[1]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_VEL_Z)) {
            total_forces[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_forces[2] * velocity[2]));
        }

        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_X)) {
            total_moment[0] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[0] * angular_velocity[0]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_Y)) {
            total_moment[1] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[1] * angular_velocity[1]));
        }
        if (central_node.IsNot(DEMFlags::FIXED_ANG_VEL_Z)) {
            total_moment[2] *= (1.0 - mGlobalDamping * GeometryFunctions::sign(total_moment[2] * angular_velocity[2]));
        }

        KRATOS_CATCH("")
    }

void PolyhedronContactElement::ClonePolyhedronDiscontinuumConstitutiveLawWithNeighbour() {
    Properties::Pointer properties_of_this_contact = mPolyhedronParticle1->GetProperties().pGetSubProperties(mPolyhedronParticle2->GetProperties().Id());
    mPolyhedronDiscontinuumConstitutiveLaw = (*properties_of_this_contact)[DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER]->CloneUnique();
	mPolyhedronDiscontinuumConstitutiveLaw->Initialize(properties_of_this_contact);
}

} // Namespace Kratos


