class Face
{
private:
	float v1, v2, v3;
public:
	Face();
	Face(int vert1, int vert2, int vert3);
	~Face();

	float getvOne();
	float getvTwo();
	float getvThree();
};