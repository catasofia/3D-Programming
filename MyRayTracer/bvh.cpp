#include "rayAccelerator.h"
#include "macros.h"

using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
			//this->n_objs = n_objs_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }


void BVH::Build(vector<Object *> &objs) {

		
			BVHNode *root = new BVHNode();

			Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			AABB world_bbox = AABB(min, max);

			for (Object* obj : objs) {
				AABB bbox = obj->GetBoundingBox();
				world_bbox.extend(bbox);
				objects.push_back(obj);
			}
			world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
			world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
			root->setAABB(world_bbox);
			nodes.push_back(root);
			build_recursive(0, objects.size(), root); // -> root node takes all the 

		}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {

		if ((right_index - left_index) <= Threshold) { // Check if the number of objects is fewer than the threshold
			node->makeLeaf(left_index, (right_index - left_index)); // Initiate current node as a leaf with primitives from objects[left_index] to objects[right_index]
		}
		else {
			// Split intersectable objects into left and right by finding a split index
			// Get largest dimension 
			AABB node_bbox = node->getAABB();
			Vector dist = node_bbox.max - node_bbox.min;

			int dimension;

			if (dist.x >= dist.y && dist.x >= dist.z) {
				dimension = 0; // x 
			}
			else if (dist.y >= dist.x && dist.y >= dist.z) {
				dimension = 1; // y 
			}
			else {
				dimension = 2; // z 
			}

			// Sort the objects on the largest dimension and find split index
			Comparator cmp;
			cmp.dimension = dimension;
			sort(objects.begin() + left_index, objects.begin() + right_index, cmp);

			// Find the split_index, where the mid point divides the primitives in a left and right side
			float mid_point = (node_bbox.max.getAxisValue(dimension) + node_bbox.min.getAxisValue(dimension)) * 0.5f ;

			// No sides empty, recalculate mid_point
			if (objects[left_index]->getCentroid().getAxisValue(dimension) > mid_point || objects[right_index - 1]->getCentroid().getAxisValue(dimension) <= mid_point) {
				mid_point = 0.0f;
				for (int i = left_index; i < right_index; i++) {
					mid_point = mid_point + objects[i]->getCentroid().getAxisValue(dimension);
				}
				mid_point = mid_point / (right_index - left_index);
			}

			int start_index = left_index;
			int end_index = right_index;
			int mid_index;
			int split_index = start_index;

			while (start_index < end_index) {
				mid_index = start_index + (end_index - start_index) * 0.5f;
				float mid_centroid_point = objects[mid_index]->getCentroid().getAxisValue(dimension);

				if (mid_centroid_point <= mid_point) {
					start_index = mid_index + 1;
					continue;
				}
				else if (mid_centroid_point > mid_point) {
					end_index = mid_index;
					continue;
				}
				break;
			}

			while (split_index < end_index) {
				if (objects[split_index]->getCentroid().getAxisValue(dimension) > mid_point) {
					break;
				}
				split_index++;
			}

			node->makeNode(nodes.end() - nodes.begin());

			// Calculate bounding boxes of left and right sides
			AABB left_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));
			AABB right_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));

			// left box until middle
			for (int i = left_index; i < split_index; i++) { 
				left_bbox.extend(objects[i]->GetBoundingBox());
			}
			// right box until end
			for (int i = split_index; i < right_index; i++) { 
				right_bbox.extend(objects[i]->GetBoundingBox());
			}

			// new child nodes
			BVHNode* left_node = new BVHNode();
			BVHNode* right_node = new BVHNode();
			left_node->setAABB(left_bbox);
			right_node->setAABB(right_bbox);

			nodes.push_back(left_node); // insert element at end
			nodes.push_back(right_node); // insert element at end

			// repeat 
			build_recursive(left_index, split_index, left_node);
			build_recursive(split_index, right_index, right_node);
		}
	}

bool BVH::Traverse(Ray& ray, Object** hit_obj, Vector& hit_point) {

			float tmp;
			float tmin = FLT_MAX;  // contains the closest primitive intersection
			bool hit = false;

			BVHNode* currentNode = nodes[0];

			// First check if we even hit the world bounding box
			if (!currentNode->getAABB().intercepts(ray, tmp)) {
				return false;
			}

			BVHNode* left_child;
			BVHNode* right_child;
			float left_dist, right_dist;
			bool left_hit, right_hit;

			while (true) {
				//If it is an inner node, run from child nodes that intercept the ray
				if (!currentNode->isLeaf()) {
					//intersection test with both child nodes
					left_child = nodes[currentNode->getIndex()];
					left_hit = left_child->getAABB().intercepts(ray, left_dist);
					if (left_child->getAABB().isInside(ray.origin)) left_dist = 0;

					right_child = nodes[currentNode->getIndex() + 1];
					right_hit = right_child->getAABB().intercepts(ray, right_dist);
					if (right_child->getAABB().isInside(ray.origin)) right_dist = 0;

					if (left_hit && right_hit) { // both nodes hit => Put the one furthest away on the stack. CurrentNode = closest node
						if (left_dist < right_dist) { //right_child on stack and currentNode left_child
							hit_stack.push(StackItem(right_child, right_dist));
							currentNode = left_child;
						}
						else { //left_child on stack and currentNode right_child
							hit_stack.push(StackItem(left_child, left_dist));
							currentNode = right_child;
						}
						continue;
					}
					else if (left_hit) { // one hit => CurrentNode = hit node
						currentNode = left_child;
						continue;
					}
					else if (right_hit) { // one hit => CurrentNode = hit node
						currentNode = right_child;
						continue;
					}
				}
				// If it is a leaf node check objects and update tmin in case of closer hit
				else {
					Object* obj;

					for (int i = currentNode->getIndex(); i < currentNode->getIndex() + currentNode->getNObjs(); i++) {
						obj = objects[i];
						if (obj->intercepts(ray, tmp) && tmp < tmin) {
							tmin = tmp;
							*hit_obj = obj;
						}
					}

					if (hit_obj != NULL) {
						hit = true;
					}
				}

				bool node_stack = false;

				// Pop stack until you find a node with t < tmin => CurrentNode = pop’d
				while (!hit_stack.empty()) {
					StackItem popped = hit_stack.top();
					hit_stack.pop();

					if (popped.t < tmin) { 
						currentNode = popped.ptr;
						node_stack = true;
						break;
					}
				}

				if (!node_stack) 
					break;
			}

			if (hit) {
				hit_point = ray.origin + ray.direction * tmin;
				return true;
			}

			return false;
}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
			float tmp;

			double length = ray.direction.length(); //distance between light and intersection point
			ray.direction.normalize();

			BVHNode* currentNode = nodes[0];

			// First check if we even hit the world bounding box
			if (!currentNode->getAABB().intercepts(ray, tmp)) {
				return false;
			}

			BVHNode* left_child;
			BVHNode* right_child;

			while (true) {
				//If it is an inner node, run from child nodes that intercept the ray
				if (!currentNode->isLeaf()) {

					float left_dist;
					left_child = nodes[currentNode->getIndex()];
					bool left_hit = left_child->getAABB().intercepts(ray, left_dist);

					float right_dist;
					right_child = nodes[currentNode->getIndex() + 1];
					bool right_hit = right_child->getAABB().intercepts(ray, right_dist);

					if (left_hit && right_hit) { // both nodes hit => Put the one furthest away on the stack. CurrentNode = closest node
						if (left_dist < right_dist) { //right_child on stack and currentNode left_child
							hit_stack.push(StackItem(right_child, right_dist));
							currentNode = left_child;
						}
						else { //left_child on stack and currentNode right_child
							hit_stack.push(StackItem(left_child, left_dist));
							currentNode = right_child;
						}
						continue;
					}
					else if (left_hit) { // one hit => CurrentNode = hit node
						currentNode = left_child;
						continue;
					}
					else if (right_hit) { // one hit => CurrentNode = hit node
						currentNode = right_child;
						continue;
					}
				}
				else {
					Object* obj;
					for (int i = currentNode->getIndex(); i < currentNode->getIndex() + currentNode->getNObjs(); i++) {
						obj = objects[i];
						if (obj->intercepts(ray, tmp) && tmp < length) {
							return true;
						}
					}
				}

				if (hit_stack.empty())
					break;
				StackItem popped = hit_stack.top();
				hit_stack.pop();
				currentNode = popped.ptr;
			}

			return false;
	}		
