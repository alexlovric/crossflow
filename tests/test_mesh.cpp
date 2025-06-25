#include "gtest/gtest.h"
#include "mesh.hpp"

TEST(MeshTest, BlockMeshCreation) {
    Mesh mesh;
    mesh.setDim(2);
    mesh.createBlockMesh2D(0, 1, 0, 1, 2, 2);
    
    EXPECT_EQ(mesh.nnp, 9) << "2x2 block mesh should have 9 nodes";
    EXPECT_EQ(mesh.nel, 4) << "2x2 block mesh should have 4 elements";
}
