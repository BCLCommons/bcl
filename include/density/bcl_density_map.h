// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_DENSITY_MAP_H_
#define BCL_DENSITY_MAP_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_nd.h"
#include "math/bcl_math_tricubic_spline.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Map
    //! @brief This class is a hashmap class
    //! @details At this point it only reads mrc files which is very difficult because there are many things you
    //! have to worry about while reading mrc files /n
    //! <A HREF="http://ami.scripps.edu/prtl_data/mrc_specification.htm"> mrc-file specifications </A>\n
    //!
    //! one of the most important things is that the bits could be swapped. There is one function,
    //! that checks this. CheckSwap( std::istream &ISTREAM) reads the first float in two directions and
    //! compares the int value. The first value of a mrcfile is the X-dimension. This is normally not so
    //! huge, so that the smaller value is considered to be the right one. Then a swap bool is that to true
    //! if the bits have to be read in the other direction
    //!
    //! @see @link example_density_map.cpp @endlink
    //! @author woetzen
    //! @date 16.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Map :
      public util::SerializableInterface
    {
    private:

    //////////
    // data //
    //////////

      // max number labels in mrc maps
      static const int s_ExtraSize       = 100;
      static const int s_MaxNumberLabels =  10;
      static const int s_LabelLength     =  80;

      linal::VectorND< int, 3> m_Index;     //!< index of map
      linal::VectorND< int, 3> m_Intervals; //!< intervals along each axis
      linal::Vector3D          m_Length;    //!< size of of unitcell in Angstrom
      linal::Vector3D          m_CellWidth; //!< size of each voxel in Angstroem
      linal::Vector3D          m_Angle;     //!< voxel angle in degrees
      linal::VectorND< int, 3> m_Axis;      //!< col, row and section axis
                                           // (columns axis: 1=X)
                                           // (rows axis: 2=Y)
                                           // (sections axis: 3=Z)
      double                     m_Minimum;  //!< Minimum intensity
      double                     m_Maximum;  //!< Maximum intensity
      double                     m_Mean;     //!< Mean    intensity
      int                        m_SpaceGroup; //!< Space group number (0 or 1)
      int                        m_NumberBytesSymmetryData; //!< number of bytes used for symmetry data 0 or 80
      char                       m_Extra[ s_ExtraSize]; //!< extra data 0 by default
      double                     m_Rmsd;     //!< rmsd of density to mean density value
      linal::Vector3D             m_Origin;   //!< origin of map
      math::Tensor< double>      m_Data;     //!< density as a tensor does also contain the dimensions of the density map

      // unused header information
      int m_MachineStamp;
      storage::Vector< std::string> m_Labels;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief default angle
      //! @return Vector3D with three default angles
      static const linal::Vector3D &GetDefaultAngle();

      //! @brief default axis
      //! @brief storage vector ND indicating indes for slow, middle and fast changing axis
      static const linal::VectorND< int, 3> &GetDefaultAxis();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Map();

      //! @brief construct from all map parameters
      //! @brief DATA the data as tensor
      Map
      (
        const math::Tensor< double> &DATA,
        const linal::VectorND< int, 3> &INDEX,
        const linal::VectorND< int, 3> &INTERVALS,
        const linal::Vector3D &LENGTH,
        const linal::Vector3D &CELLWIDTH,
        const linal::Vector3D &ANGLE,
        const linal::VectorND< int, 3> &AXIS,
        const linal::Vector3D &ORIGIN
      );

      //! @brief copy constructor
      //! @param DENSITY_MAP map to copy from
      Map( const Map &DENSITY_MAP);

      //! copy constructor
      Map *Clone() const
      {
        return new Map( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "DensityMap");
        return s_name;
      }

      //! return number of Elements
      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! return Index
      const linal::VectorND< int, 3> &GetIndex() const
      {
        return m_Index;
      }

      //! return Intervals
      const linal::VectorND< int, 3> &GetIntervals() const
      {
        return m_Intervals;
      }

      //! return Dimensions
      storage::VectorND< 3, size_t> GetDimensions() const
      {
        return storage::VectorND< 3, size_t>( m_Data.GetNumberCols(), m_Data.GetNumberRows(), m_Data.NumberLayers());
      }

      //! return Width of unitcell
      const linal::Vector3D &GetUnitCellLength() const
      {
        return m_Length;
      }

      //! return Width of each cell
      const linal::Vector3D &GetCellWidth() const
      {
        return m_CellWidth;
      }

      //! return Angles of unit cell
      const linal::Vector3D &GetAngles() const
      {
        return m_Angle;
      }

      //! return Angles of unit cell
      const linal::VectorND< int, 3> &GetAxis() const
      {
        return m_Axis;
      }

      //! return Origin for Transformations
      const linal::Vector3D &GetOrigin() const
      {
        return m_Origin;
      }

      //! return reference to Value (x, y, z)
      double const &operator()( const size_t COL, const size_t ROW, const size_t LAYER) const
      {
        return m_Data( LAYER, ROW, COL);
      }

      //! return reference to changeable Value (x, y, z)
      double &operator()( const size_t COL, const size_t ROW, const size_t LAYER)
      {
        return m_Data( LAYER, ROW, COL);
      }

      //! return reference to changeable Value (x, y, z)
      double &operator()( const storage::VectorND< 3, size_t> &INDEX)
      {
        return m_Data( INDEX.Third(), INDEX.Second(), INDEX.First());
      }

      //! return reference to Value (x, y, z)
      double const &operator()( const storage::VectorND< 3, size_t> &INDEX) const
      {
        return m_Data( INDEX.Third(), INDEX.Second(), INDEX.First());
      }

      //! return reference to changeable Value (x, y, z)
      double &operator()( const linal::VectorND< size_t, 3> &INDEX)
      {
        return m_Data( INDEX( 2), INDEX( 1), INDEX( 0));
      }

      //! return reference to changeable Value (x, y, z)
      const double &operator()( const linal::VectorND< size_t, 3> &INDEX) const
      {
        return m_Data( INDEX( 2), INDEX( 1), INDEX( 0));
      }

      //! return minimal value of Map
      double GetMinimum() const
      {
        return m_Minimum;
      }

      //! return maximal value of Map
      double GetMaximum() const
      {
        return m_Maximum;
      }

      //! return mean value of Map
      double GetMean() const
      {
        return m_Mean;
      }

      //! return Rmsd to Mean value in density map
      double GetRmsd() const
      {
        return m_Rmsd;
      }

      //! @brief calculate the volume of one voxel
      //! @return the volume of a voxel V = cellwidth x * y * z
      double GetVoxelVolume() const
      {
        return m_CellWidth.X() * m_CellWidth.Y() * m_CellWidth.Z();
      }

      //! return const reference to the Tensor density map
      const math::Tensor< double> &GetData() const
      {
        return m_Data;
      }

      //! return reference to changeable Tensor of the density map
      math::Tensor< double> &GetData()
      {
        return m_Data;
      }

      //! @brief volume of density map
      //! @return the volume determined by size*voxelvolume
      double GetVolume() const
      {
        return m_Data.GetSize() * GetVoxelVolume();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the min max and mean from the intensities stored in the data
      void CalculateMinMaxMeanRmsd();

      //! compute histogram over all densities in NR_BINS bins
      //! @param NR_BINS number of bins in histogram
      //! @return Histogram
      math::Histogram Histogram( const size_t NR_BINS) const;

      //! cut out a part of the density map
      //! @param POSCOL   column or x position
      //! @param POSROW   row    pr y position
      //! @param POSLAYER layer  or z position
      //! @param NCOL     number of columns
      //! @param NROW     number of rows
      //! @param NLAYER   number of layers
      //! @return sub density map from given position of given size
      Map SubMap
      (
        const size_t POSCOL, const size_t POSROW, const size_t POSLAYER,
        const size_t NCOL, const size_t NROW, const size_t NLAYER
      ) const;

      //! get the common sub tensor from this and argument density map
      //! @param DENSITY_MAP the second density map to find common subtensor
      //! @return two tensors, first tensor of this map, second tensot of argument
      storage::VectorND< 2, math::Tensor< double> > CommonSubTensor( const Map &DENSITY_MAP) const;

      //! compute a point cloud of NUMBER_OF_POINTS with an at least distance of FEATURE_DISTANCE in Angstrom using
      //! the highest derivatives between voxels in the map and the own intensity as the user can weight
      coord::PointCloud CalculatePointCloud
      (
        const size_t NUMBER_OF_POINTS,
        const double FEATURE_DISTANCE,
        const double RATIO_INTENSITY_GRADIENT
      ) const;

      //! @brief add noise
      //! @param RNG an random number Generator object
      //! @param MEAN the mean of the Gaussian distribution
      //! @param STANDARD_DEVIATION standard deviation of the Gaussian distribution
      //! @return the cross correlation coefficient map (current map as simulated) vs new map with noise
      double AddNoise( const random::DistributionInterface &RNG, const double MEAN, const double STANDARD_DEVIATION);

      //! @brief normalize
      //! transforms intensities, so that rmsd is 1
      void Normalize();

      //! calculate Correlation factor to a density map
      double Correlation( const Map &DENSITY_MAP) const;

      //! calculate the cross correlation factor to a simulated density map, only using voxels that are above the given
      //! CONTOUR_LEVEL in the simulated density map
      //! @param SIMULATED_DENSITY_MAP
      //! @param CONTOUR_LEVEL
      //! @return cross correlation coefficient
      double CrossCorrelationCoefficient
      (
        const Map &SIMULATED_DENSITY_MAP, const double CONTOUR_LEVEL
      ) const;

      //! calculate the standard deviation to a smaller density map
      double StandardDeviation( const Map &DENSITY_MAP) const;

      //! @brief determine edges using a Sobel operator
      //! @sa http://en.wikipedia.org/wiki/Edge_detection
      Map &EdgeDetectionSobel
      (
        const double WEIGTH_FACE = double( 6),
        const double WEIGHT_EDGE = double( 3),
        const double WEIGHT_VERTEX = double( 1)
      );

      //! @brief balance the intensities
      //! @param LENGTH of cube vertex in which intensities are balanced
      Map &BalanceIntensities( const double LENGTH);

      //! @brief convert Map into Spline
      math::TricubicSpline ConvertDensityToSpline() const;

      //! @brief orthogonalize map if angles aren't 90,90,90
      Map OrthogonalizeMap() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param MAP_RHS right hand side density map to be copied
      //! @return Map the changed map
      Map &operator=( const Map &MAP_RHS);

    //////////////////////
    // input and output //
    //////////////////////

      //! write Map to std::ostream using the given util::Format
      std::ostream &WriteHeader( std::ostream &OSTREAM) const;

      //! read Map from mrc file EXTENDED_HEADER is the number of bytes additional to a normal 1024bytes header
      std::istream &ReadMRC( std::istream &ISTREAM, const size_t EXTENDED_HEADER = 0);

      //! write Map to mrc file
      std::ostream &WriteMRC( std::ostream &OSTREAM) const;

    protected:

      //! write Map to std::ostream using standard util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read Map from io::IFStream
      std::istream &Read( std::istream &ISTREAM);

    private:

      //!checks stream of mrc file whether the bytes are swapped
      static bool CheckSwap( std::istream &ISTREAM);

      //!read int from istream and swap if necessary
      static int ReadInt( std::istream &ISTREAM, const bool SWAP);

      //!read float from istream and swap if necessary
      static float ReadFloat( std::istream &ISTREAM, const bool SWAP);

      //!order values according to given Axis order
      template< typename t_DataType>
      static linal::VectorND< t_DataType, 3>
      OrderValues( const linal::VectorND< int, 3> &AXIS, const linal::VectorConstInterface< t_DataType> &VALUES)
      {
        linal::VectorND< int, 3> index_xyz;
        index_xyz( AXIS( 0) - 1) = 0;
        index_xyz( AXIS( 1) - 1) = 1;
        index_xyz( AXIS( 2) - 1) = 2;

        linal::VectorND< t_DataType, 3> new_values;
        new_values( 0) = VALUES( index_xyz( 0));
        new_values( 1) = VALUES( index_xyz( 1));
        new_values( 2) = VALUES( index_xyz( 2));

        return new_values;
      }

      //!order values according to given Axis order
      static linal::Vector3D OrderValues( const linal::VectorND< int, 3> &AXIS, const linal::Vector3D &VALUES);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class Voxel
      //! @brief helper struct for conversion of density map to a pointcloud
      //! @author woetzen
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Voxel
      {
        size_t m_Layer;     //!< Pos Layer
        size_t m_Row;       //!< Pos Row
        size_t m_Col;       //!< Pos Col
        double m_Intensity; //!< Intensity
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class VoxelDensityGreater
      //! @brief compare voxels by there intensity
      //! @author woetzen
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct VoxelDensityGreater :
        public std::binary_function< Voxel, Voxel, bool>
      {
        bool operator()( const Voxel &VOXEL_LHS, const Voxel &VOXEL_RHS) const
        {
          return VOXEL_LHS.m_Intensity > VOXEL_RHS.m_Intensity;
        }
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class VoxelOverlap
      //! @brief compare voxels if they are within each others radius
      //! @author woetzen
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct VoxelOverlap :
        public std::binary_function< Voxel, Voxel, bool>
      {
        const double m_SquareDistance;
        const size_t m_ThreshX;
        const size_t m_ThreshY;
        const size_t m_ThreshZ;
        const linal::Vector3D m_CellWidth;

        //! @brief construct from voxel property and the distance for overlap
        //! @param FEATURE_DISTANCE - if two voxels are found to be within that distance, they are considered overlapping
        //! @param CELL_WIDTH the side length of the supplied voxels
        VoxelOverlap
        (
          const double FEATURE_DISTANCE,
          const linal::Vector3D &CELL_WIDTH
        ) :
          m_SquareDistance( FEATURE_DISTANCE * FEATURE_DISTANCE),
          m_ThreshX( size_t( std::ceil( FEATURE_DISTANCE / CELL_WIDTH.X()))),
          m_ThreshY( size_t( std::ceil( FEATURE_DISTANCE / CELL_WIDTH.Y()))),
          m_ThreshZ( size_t( std::ceil( FEATURE_DISTANCE / CELL_WIDTH.Z()))),
          m_CellWidth( CELL_WIDTH)
        {
        }

        //! @brief operator that checks if two voxels overlap given a minimal distance and a cell width
        //! @param VOXEL_LHS left hand side voxel
        //! @param VOXEL_RHS right hand side svoxel
        //! @return true if the two voxel overlap considering the middle point of the voxel and the and the m_MinimalDistance
        bool operator()( const Voxel &VOXEL_LHS, const Voxel &VOXEL_RHS) const
        {
          const size_t delta_x( math::Absolute( int( VOXEL_LHS.m_Col)   - int( VOXEL_RHS.m_Col)));
          const size_t delta_y( math::Absolute( int( VOXEL_LHS.m_Row)   - int( VOXEL_RHS.m_Row)));
          const size_t delta_z( math::Absolute( int( VOXEL_LHS.m_Layer) - int( VOXEL_RHS.m_Layer)));

          if( delta_x <= m_ThreshX || delta_y <= m_ThreshY || delta_z <= m_ThreshZ)
          {
            if
            (
              math::Sqr( m_CellWidth.X() * delta_x) +
              math::Sqr( m_CellWidth.Y() * delta_y) +
              math::Sqr( m_CellWidth.Z() * delta_z)
              < m_SquareDistance
            )
            {
              return true;
            }
          }

          // does not overlap
          return false;
        }

      }; // struct VoxelOverlap

    }; // class Map

  } // namespace density
} // namespace bcl

#endif //BCL_DENSITY_MAP_H_
