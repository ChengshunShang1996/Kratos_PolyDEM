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

    mIsFaceParallel = false;

	if (GJK()) {
        
		auto& central_node_1 = mPolyhedronParticle1->GetGeometry()[0];
		auto& central_node_2 = mPolyhedronParticle2->GetGeometry()[0];
		Vector3 coll1Pos = {central_node_1[0], central_node_1[1], central_node_1[2]};
		Vector3 coll2Pos = {central_node_2[0], central_node_2[1], central_node_2[2]};

		Vector3 contact_force(0.0, 0.0, 0.0);
		Vector3 contact_moment_1(0.0, 0.0, 0.0);
		Vector3 contact_moment_2(0.0, 0.0, 0.0);

		//double kn = 100000.0;
		//contact_force = mOverlapVector * kn;

		ClonePolyhedronDiscontinuumConstitutiveLawWithNeighbour();
    	mPolyhedronDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, mPolyhedronParticle1, mPolyhedronParticle2, 
																mOverlapVector, mContactPoint, contact_force, mTangentialElasticContactForce);

		if (r_process_info[ROTATION_OPTION]){
			Vector3 torque_arm_1 = mContactPoint - coll1Pos;
			Vector3 torque_arm_2 = mContactPoint - coll2Pos;
			contact_moment_1 = Vector3::Cross(torque_arm_1, contact_force);
			contact_moment_2 = Vector3::Cross(torque_arm_2, -contact_force);
			//KRATOS_INFO("torque_arm_1") << torque_arm_1 << std::endl;
			//KRATOS_INFO("torque_arm_2") << torque_arm_2 << std::endl;
			//KRATOS_INFO("Contact Moment 1") << contact_moment_1 << std::endl;
			//KRATOS_INFO("Contact Moment 2") << contact_moment_2 << std::endl;
			const array_1d<double, 3>& velocity_1 = central_node_1.FastGetSolutionStepValue(VELOCITY);
        	const array_1d<double, 3>& velocity_2 = central_node_2.FastGetSolutionStepValue(VELOCITY);
			Vector3 velocity_1_vec(velocity_1[0], velocity_1[1], velocity_1[2]);
			Vector3 velocity_2_vec(velocity_2[0], velocity_2[1], velocity_2[2]);
			double velocity_1_vec_modulus = velocity_1_vec.Length();
			double velocity_2_vec_modulus = velocity_2_vec.Length();
			if (velocity_1_vec_modulus > 1e-3 && mIsFaceParallel == true){
				central_node_1.FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] *= 0.5;
				central_node_1.FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] *= 0.5;
				central_node_1.FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] *= 0.5;
			}
			if (velocity_2_vec_modulus > 1e-3 && mIsFaceParallel == true){
				central_node_2.FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] *= 0.5;
				central_node_2.FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] *= 0.5;
				central_node_2.FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] *= 0.5;
			}
		}

		array_1d<double,3>& contact_forces_1 = central_node_1.FastGetSolutionStepValue(CONTACT_FORCES);
		array_1d<double,3>& total_forces_1 = central_node_1.FastGetSolutionStepValue(TOTAL_FORCES);
		array_1d<double,3>& total_moment_1 = central_node_1.FastGetSolutionStepValue(PARTICLE_MOMENT);

		//TODO: can be improved by calculating the weighted contact forces
		contact_forces_1[0] = contact_force[0];
		contact_forces_1[1] = contact_force[1];
		contact_forces_1[2] = contact_force[2];
		
		total_forces_1[0] += contact_force[0];
		total_forces_1[1] += contact_force[1];
		total_forces_1[2] += contact_force[2];

		total_moment_1[0] += contact_moment_1[0];
		total_moment_1[1] += contact_moment_1[1];
		total_moment_1[2] += contact_moment_1[2];

		array_1d<double,3>& contact_forces_2 = central_node_2.FastGetSolutionStepValue(CONTACT_FORCES);
		array_1d<double,3>& total_forces_2 = central_node_2.FastGetSolutionStepValue(TOTAL_FORCES);
		array_1d<double,3>& total_moment_2 = central_node_2.FastGetSolutionStepValue(PARTICLE_MOMENT);

		contact_forces_2[0] = -contact_force[0];
		contact_forces_2[1] = -contact_force[1];
		contact_forces_2[2] = -contact_force[2];
		
		total_forces_2[0] -= contact_force[0];
		total_forces_2[1] -= contact_force[1];
		total_forces_2[2] -= contact_force[2];

		total_moment_2[0] += contact_moment_2[0];
		total_moment_2[1] += contact_moment_2[1];
		total_moment_2[2] += contact_moment_2[2];

		//KRATOS_INFO("Contact Force") << contact_force << std::endl;
		//KRATOS_INFO("TOTAL_FORCES_1") << total_forces_1 << std::endl;

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
			//KRATOS_WATCH(mOverlapVector);
			//Vector3 normal = (mContactPoint1 - mContactPoint2).Normalised();
            //Vector3 contact_normal = faces[closest_face][3]
            //Vector3 contact_point_minkowski = contact_normal * Vector3::Dot(p.p, search_dir)

			// Identify the corresponding geometric features on the two polyhedra
			bool poly1_find_face = true;
			bool poly2_find_face = true;
			std::vector<Vector3> faceVertices1 = mPolyhedronParticle1->GetIntersectingFaceVertices(mContactPoint1, mOverlapVector, poly1_find_face);
			std::vector<Vector3> faceVertices2 = mPolyhedronParticle2->GetIntersectingFaceVertices(mContactPoint2, -mOverlapVector, poly2_find_face);
			
			//if (faceVertices1.size() == 0 || faceVertices2.size() == 0) {
			//	mContactPoint = (mContactPoint1 + mContactPoint2)/2;
			//	return;
			//}
			
			Vector3 normal1 = CalculateFaceNormal(faceVertices1);
			Vector3 normal2 = CalculateFaceNormal(faceVertices2);

			// Check if it is a face-face contact
			if (AreNormalsParallel(normal1, normal2) && poly1_find_face && poly2_find_face) {
				// Calculate the contact area for face-face contact
				std::vector<Vector3> intersectionVertices = CalculateIntersection(faceVertices1, faceVertices2, normal1);

				if (!intersectionVertices.empty()){
					// Determine the contact point
					mContactPoint = CalculateCentroid(intersectionVertices);

					Vector3 TempOverlapVector = mOverlapVector.Normalised();

					if (Vector3::Dot(normal1, TempOverlapVector) < 0){
						normal1	= -normal1;
					} 

					mOverlapVector = normal1 * mOverlapVector.Length();

					mIsFaceParallel = true;

				} else {
					mContactPoint = (mContactPoint1 + mContactPoint2)/2;
				}

				return;
			} else {

				mContactPoint = (mContactPoint1 + mContactPoint2)/2;

				return;
			}
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
	KRATOS_ERROR << "EPA did not converge" << std::endl;
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

	// Identify the corresponding geometric features on the two polyhedra
	bool poly1_find_face = true;
	bool poly2_find_face = true;
	std::vector<Vector3> faceVertices1 = mPolyhedronParticle1->GetIntersectingFaceVertices(mContactPoint1, mOverlapVector, poly1_find_face);
	std::vector<Vector3> faceVertices2 = mPolyhedronParticle2->GetIntersectingFaceVertices(mContactPoint2, -mOverlapVector, poly2_find_face);
	
	//if (faceVertices1.size() == 0 || faceVertices2.size() == 0) {
	//	mContactPoint = (mContactPoint1 + mContactPoint2)/2;
	//	return;
	//}
	
	Vector3 normal1 = CalculateFaceNormal(faceVertices1);
	Vector3 normal2 = CalculateFaceNormal(faceVertices2);

	// Check if it is a face-face contact
	if (AreNormalsParallel(normal1, normal2) && poly1_find_face && poly2_find_face) {
		// Calculate the contact area for face-face contact
		std::vector<Vector3> intersectionVertices = CalculateIntersection(faceVertices1, faceVertices2, normal1);

		if (!intersectionVertices.empty()){
			// Determine the contact point
			mContactPoint = CalculateCentroid(intersectionVertices);

			Vector3 TempOverlapVector = mOverlapVector.Normalised();

			if (Vector3::Dot(normal1, TempOverlapVector) < 0){
				normal1	= -normal1;
			} 

			mOverlapVector = normal1 * mOverlapVector.Length();

			mIsFaceParallel = true;

		} else {
			mContactPoint = (mContactPoint1 + mContactPoint2)/2;
		}
		
		return;

	} else {

		mContactPoint = (mContactPoint1 + mContactPoint2)/2;

		return;
	}

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
	if (std::abs(denom) < std::numeric_limits<double>::epsilon()) {
        u = v = w = 0.0; 
        return;
    }
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;

    KRATOS_CATCH( "" )
}

bool PolyhedronContactElement::AreNormalsParallel(const Vector3& normal1, const Vector3& normal2)
{
    // Normalize the normals
    Vector3 norm1 = normal1.Normalised();
    Vector3 norm2 = normal2.Normalised();

    // Calculate the dot product
    double dotProduct = Vector3::Dot(norm1, norm2);

    // Check if the dot product is close to 1 or -1
    const double epsilon = 1e-6; // Tolerance for floating-point comparison
    return std::abs(std::abs(dotProduct) - 1.0) < epsilon;
}

Vector3 PolyhedronContactElement::CalculateFaceNormal(const std::vector<Vector3>& vertices)
{
	// Ensure there are at least 3 vertices to define a plane
    if (vertices.size() < 3) {
        KRATOS_INFO("ERROR") << "A face must have at least 3 vertices to calculate a normal." << std::endl;
    }
    // Calculate the normal using the cross product of two edges of the face
    Vector3 edge1 = vertices[1] - vertices[0];
    Vector3 edge2 = vertices[2] - vertices[0];
    Vector3 normal = Vector3::Cross(edge1, edge2);

    // Normalize the normal vector
    normal.Normalise();

    return normal;
}

std::vector<Vector3> PolyhedronContactElement::CalculateIntersection(const std::vector<Vector3>& faceVertices1, const std::vector<Vector3>& faceVertices2, const Vector3& normal)
{
    // Project the vertices of both faces onto a common plane defined by the normal
    std::vector<Vector2> projectedFace1 = ProjectToPlane(faceVertices1, normal);
    std::vector<Vector2> projectedFace2 = ProjectToPlane(faceVertices2, normal);

	Vector3 u;
	if (std::abs(normal[0]) > std::abs(normal[1])) {
		u = Vector3(-normal[2], 0, normal[0]).Normalised();
	} else {
		u = Vector3(0, -normal[2], normal[1]).Normalised();
	}
    Vector3 v = Vector3::Cross(normal, u);

	Vector3 projectedPoint3d = u * projectedFace1[0].x + v * projectedFace1[0].y;
    Vector3 relativeVector = faceVertices1[0] - projectedPoint3d;

    // Calculate the intersection of the two projected polygons
    std::vector<Vector2> intersection2D = CalculatePolygonIntersection(projectedFace1, projectedFace2);

	// Reproject the intersection vertices back to 3D space
	std::vector<Vector3> intersection3D = ReprojectTo3D(intersection2D, normal, relativeVector);

    return intersection3D;
}

std::vector<Vector2> PolyhedronContactElement::ProjectToPlane(const std::vector<Vector3>& vertices, const Vector3& normal)
{
    // Choose a basis for the plane
    Vector3 u;
	if (std::abs(normal[0]) > std::abs(normal[1])) {
		u = Vector3(-normal[2], 0, normal[0]).Normalised();
	} else {
		u = Vector3(0, -normal[2], normal[1]).Normalised();
	}
    Vector3 v = Vector3::Cross(normal, u);

    // Project each vertex onto the plane
    std::vector<Vector2> projectedVertices;
    for (const auto& vertex : vertices) {
        double x = Vector3::Dot(vertex, u);
        double y = Vector3::Dot(vertex, v);
        projectedVertices.push_back(Vector2(x, y));
    }

    return projectedVertices;
}

std::vector<Vector3> PolyhedronContactElement::ReprojectTo3D(const std::vector<Vector2>& vertices2D, const Vector3& normal, const Vector3& reletiveVector)
{
	Vector3 u;
	if (std::abs(normal[0]) > std::abs(normal[1])) {
		u = Vector3(-normal[2], 0, normal[0]).Normalised();
	} else {
		u = Vector3(0, -normal[2], normal[1]).Normalised();
	}
    Vector3 v = Vector3::Cross(normal, u);

    // Reproject each 2D vertex back to 3D space
    std::vector<Vector3> vertices3D;
    for (const auto& vertex2D : vertices2D) {
        Vector3 vertex3D = reletiveVector + u * vertex2D.x + v * vertex2D.y;
        vertices3D.push_back(vertex3D);
    }

    return vertices3D;
}

// Helper function to calculate the signed area of a polygon
double PolyhedronContactElement::calculateSignedArea(const std::vector<Vector2>& polygon) {
	double area = 0.0;
	size_t n = polygon.size();
	for (size_t i = 0; i < n; ++i) {
		const Vector2& p1 = polygon[i];
		const Vector2& p2 = polygon[(i + 1) % n];
		area += (p2.x - p1.x) * (p2.y + p1.y);
	}
	return area;
}

// Check if a polygon is clockwise
bool PolyhedronContactElement::isClockwise(const std::vector<Vector2>& polygon) {
	return calculateSignedArea(polygon) > 0;
}

// Ensure that a polygon is in clockwise order
void PolyhedronContactElement::ensureClockwise(std::vector<Vector2>& polygon) {
	if (!isClockwise(polygon)) {
		reverse(polygon.begin(), polygon.end());
	}
}

std::vector<Vector2> PolyhedronContactElement::CalculatePolygonIntersection(std::vector<Vector2>& subjectPolygon, std::vector<Vector2>& clippingPolygon) {
	
	ensureClockwise(subjectPolygon);
    ensureClockwise(clippingPolygon);
	
	std::vector<Vector2> finalPolygon = subjectPolygon;

	for (size_t i = 0; i < clippingPolygon.size(); ++i) {
		std::vector<Vector2> nextPolygon = finalPolygon;
		finalPolygon.clear();

		Vector2 c_edge_start = clippingPolygon[i == 0 ? clippingPolygon.size() - 1 : i - 1];
		Vector2 c_edge_end = clippingPolygon[i];

		for (size_t j = 0; j < nextPolygon.size(); ++j) {
			Vector2 s_edge_start = nextPolygon[j == 0 ? nextPolygon.size() - 1 : j - 1];
			Vector2 s_edge_end = nextPolygon[j];

			if (IsInside(c_edge_start, c_edge_end, s_edge_end)) {
				if (!IsInside(c_edge_start, c_edge_end, s_edge_start)) {
					finalPolygon.push_back(ComputeIntersection(s_edge_start, s_edge_end, c_edge_start, c_edge_end));
				}
				finalPolygon.push_back(s_edge_end);
			} else if (IsInside(c_edge_start, c_edge_end, s_edge_start)) {
				finalPolygon.push_back(ComputeIntersection(s_edge_start, s_edge_end, c_edge_start, c_edge_end));
			}
		}
	}

	return finalPolygon;
}

bool PolyhedronContactElement::IsInside(const Vector2& edgeStart, const Vector2& edgeEnd, const Vector2& point)
{
    double R = (edgeEnd.x - edgeStart.x) * (point.y - edgeStart.y) > (edgeEnd.y - edgeStart.y) * (point.x - edgeStart.x);
    return R <= 0;
}

Vector2 PolyhedronContactElement::ComputeIntersection(const Vector2& p1, const Vector2& p2, const Vector2& p3, const Vector2& p4) {
	
	double x, y;

	// If first line is vertical
	if (std::abs(p2.x - p1.x) < 1e-9) {
		x = p1.x;
		double m2 = (p4.y - p3.y) / (p4.x - p3.x);
		double b2 = p3.y - m2 * p3.x;
		y = m2 * x + b2;
	}
	// If second line is vertical
	else if (std::abs(p4.x - p3.x) < 1e-9) {
		x = p3.x;
		double m1 = (p2.y - p1.y) / (p2.x - p1.x);
		double b1 = p1.y - m1 * p1.x;
		y = m1 * x + b1;
	}
	// If neither line is vertical
	else {
		double m1 = (p2.y - p1.y) / (p2.x - p1.x);
		double b1 = p1.y - m1 * p1.x;
		double m2 = (p4.y - p3.y) / (p4.x - p3.x);
		double b2 = p3.y - m2 * p3.x;

		x = (b2 - b1) / (m1 - m2);
		y = m1 * x + b1;
	}

	return Vector2(x, y);
}

Vector3 PolyhedronContactElement::CalculateCentroid(const std::vector<Vector3>& vertices)
{
    Vector3 centroid(0.0, 0.0, 0.0);
    for (const auto& vertex : vertices) {
        centroid += vertex;
    }
    
	centroid /= vertices.size();
	
	return centroid;
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


