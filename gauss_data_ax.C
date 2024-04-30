//DYNELA - PS

    4, // Number of integration points of the Element
    {
        // Integration point 1
        {
            Vec3D(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 0.0),                            // Coordinates of the integration point
            1,                                                                         // Weight associated to the integration point
            Vector(4, (2.0 + sqrt(3.0)) / 6, 1.0 / 6, (2.0 - sqrt(3.0)) / 6, 1.0 / 6), // Shape functions at the integration point
            Matrix(4, 2,                                                               // Derivatives of the shape functions at the integration point
                   (-3.0 - sqrt(3.0)) / 12, (-3.0 - sqrt(3.0)) / 12,                   //
                   (+3.0 + sqrt(3.0)) / 12, (-3.0 + sqrt(3.0)) / 12,                   //
                   (+3.0 - sqrt(3.0)) / 12, (+3.0 - sqrt(3.0)) / 12,                   //
                   (-3.0 + sqrt(3.0)) / 12, (+3.0 + sqrt(3.0)) / 12)                   //
        },

//DYNELA - AXIL 
        // Integration point 1
        {
            Vec3D(-1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 0.0),                            // Vec3D coords (Coordinates of the integration point)
            1,                                                                         // double weight (Weight associated to the integration point)
            Vector(4, (2.0 + sqrt(3.0)) / 6, 1.0 / 6, (2.0 - sqrt(3.0)) / 6, 1.0 / 6), // Vector shapeFunction (Shape functions at the integration point)
            Matrix(4, 2,                                                               // Matrix derShapeFunction (Derivatives of the shape functions at the integration point)
                   (-3.0 - sqrt(3.0)) / 12, (-3.0 - sqrt(3.0)) / 12,                   //
                   (+3.0 + sqrt(3.0)) / 12, (-3.0 + sqrt(3.0)) / 12,                   //
                   (+3.0 - sqrt(3.0)) / 12, (+3.0 - sqrt(3.0)) / 12,                   //
                   (-3.0 + sqrt(3.0)) / 12, (+3.0 + sqrt(3.0)) / 12)                   //
        },
        // Integration point 2
        {
            Vec3D(+1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 0.0),                            // Vec3D coords (Coordinates of the integration point)
            1,                                                                         // double weight (Weight associated to the integration point)
            Vector(4, 1.0 / 6, (2.0 + sqrt(3.0)) / 6, 1.0 / 6, (2.0 - sqrt(3.0)) / 6), // Vector shapeFunction (Shape functions at the integration point)
            Matrix(4, 2,                                                               // Matrix derShapeFunction (Derivatives of the shape functions at the integration point)
                   (-3.0 - sqrt(3.0)) / 12, (-3.0 + sqrt(3.0)) / 12,                   //
                   (+3.0 + sqrt(3.0)) / 12, (-3.0 - sqrt(3.0)) / 12,                   //
                   (+3.0 - sqrt(3.0)) / 12, (+3.0 + sqrt(3.0)) / 12,                   //
                   (-3.0 + sqrt(3.0)) / 12, (+3.0 - sqrt(3.0)) / 12)                   //
        },
        // Integration point 3
        {
            Vec3D(+1.0 / sqrt(3.0), +1.0 / sqrt(3.0), 0.0),                            // Vec3D coords (Coordinates of the integration point)
            1,                                                                         // double weight (Weight associated to the integration point)
            Vector(4, (2.0 - sqrt(3.0)) / 6, 1.0 / 6, (2.0 + sqrt(3.0)) / 6, 1.0 / 6), // Vector shapeFunction (Shape functions at the integration point)
            Matrix(4, 2,                                                               // Matrix derShapeFunction (Derivatives of the shape functions at the integration point)
                   (-3.0 + sqrt(3.0)) / 12, (-3.0 + sqrt(3.0)) / 12,                   //
                   (+3.0 - sqrt(3.0)) / 12, (-3.0 - sqrt(3.0)) / 12,                   //
                   (+3.0 + sqrt(3.0)) / 12, (+3.0 + sqrt(3.0)) / 12,                   //
                   (-3.0 - sqrt(3.0)) / 12, (+3.0 - sqrt(3.0)) / 12)                   //
        },
        // Integration point 4
        {
            Vec3D(-1.0 / sqrt(3.0), +1.0 / sqrt(3.0), 0.0),                            // Vec3D coords (Coordinates of the integration point)
            1,                                                                         // double weight (Weight associated to the integration point)
            Vector(4, 1.0 / 6, (2.0 - sqrt(3.0)) / 6, 1.0 / 6, (2.0 + sqrt(3.0)) / 6), // Vector shapeFunction (Shape functions at the integration point)
            Matrix(4, 2,                                                               // Matrix derShapeFunction (Derivatives of the shape functions at the integration point)
                   (-3.0 + sqrt(3.0)) / 12, (-3.0 - sqrt(3.0)) / 12,                   //
                   (+3.0 - sqrt(3.0)) / 12, (-3.0 + sqrt(3.0)) / 12,                   //
                   (+3.0 + sqrt(3.0)) / 12, (+3.0 - sqrt(3.0)) / 12,                   //
                   (-3.0 - sqrt(3.0)) / 12, (+3.0 + sqrt(3.0)) / 12)                   //
        }
        //
    },