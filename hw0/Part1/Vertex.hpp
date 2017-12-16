class Vertex
{
private:
	float x, y, z;
public:
	Vertex();
	Vertex(float x1, float x2, float x3);
	~Vertex();

	float getX();
	float getY();
	float getZ();
};