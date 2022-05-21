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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "density/bcl_density_map.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_cloud.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Map::s_Instance
    (
      GetObjectInstances().AddInstance( new Map())
    );

    //! @brief default angle
    //! @return Vector3D with three default angles
    const linal::Vector3D &Map::GetDefaultAngle()
    {
      // default is 90 degrees
      static const linal::Vector3D s_angle( 90.0, 90.0, 90.0);
      return s_angle;
    }

    //! @brief default axis
    //! @brief storage vector ND indicating indes for slow, middle and fast changing axis
    const linal::VectorND< int, 3> &Map::GetDefaultAxis()
    {
      // default is first index for slow changing index, 3 for fast chaning index
      static const linal::VectorND< int, 3> s_axis( 1, 2, 3);
      return s_axis;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Map::Map() :
      m_Index(),
      m_Intervals(),
      m_Length(),
      m_CellWidth(),
      m_Angle( GetDefaultAngle()),
      m_Axis( GetDefaultAxis()),
      m_Minimum(),
      m_Maximum(),
      m_Mean(),
      m_SpaceGroup( 0),
      m_NumberBytesSymmetryData( 0),
      m_Rmsd(),
      m_Origin(),
      m_Data(),
      m_MachineStamp( 0)
    {
      std::fill_n( m_Extra, s_ExtraSize, '0');
    }

    //! @brief construct from all map parameters
    //! @brief DATA the data as tensor
    Map::Map
    (
      const math::Tensor< double> &DATA,
      const linal::VectorND< int, 3> &INDEX,
      const linal::VectorND< int, 3> &INTERVALS,
      const linal::Vector3D &LENGTH,
      const linal::Vector3D &CELLWIDTH,
      const linal::Vector3D &ANGLE,
      const linal::VectorND< int, 3> &AXIS,
      const linal::Vector3D &ORIGIN
    ) :
      m_Index( INDEX),
      m_Intervals( INTERVALS),
      m_Length( LENGTH),
      m_CellWidth( CELLWIDTH),
      m_Angle( ANGLE),
      m_Axis( AXIS),
      m_Minimum(),
      m_Maximum(),
      m_Mean(),
      m_SpaceGroup( 0),
      m_NumberBytesSymmetryData( 0),
      m_Rmsd(),
      m_Origin( ORIGIN),
      m_Data( DATA),
      m_MachineStamp( 0),
      m_Labels()
    {
      std::fill_n( m_Extra, s_ExtraSize, '0');
      CalculateMinMaxMeanRmsd();
    }

    //! @brief copy constructor
    //! @param DENSITY_MAP map to copy from
    Map::Map( const Map &DENSITY_MAP) :
      m_Index( DENSITY_MAP.m_Index),
      m_Intervals( DENSITY_MAP.m_Intervals),
      m_Length( DENSITY_MAP.m_Length),
      m_CellWidth( DENSITY_MAP.m_CellWidth),
      m_Angle( DENSITY_MAP.m_Angle),
      m_Axis( DENSITY_MAP.m_Axis),
      m_Minimum( DENSITY_MAP.m_Minimum),
      m_Maximum( DENSITY_MAP.m_Maximum),
      m_Mean( DENSITY_MAP.m_Mean),
      m_SpaceGroup( DENSITY_MAP.m_SpaceGroup),
      m_NumberBytesSymmetryData( DENSITY_MAP.m_NumberBytesSymmetryData),
      m_Rmsd( DENSITY_MAP.m_Rmsd),
      m_Origin( DENSITY_MAP.m_Origin),
      m_Data( DENSITY_MAP.m_Data),
      m_MachineStamp( DENSITY_MAP.m_MachineStamp),
      m_Labels( DENSITY_MAP.m_Labels)
    {
      std::copy( DENSITY_MAP.m_Extra, DENSITY_MAP.m_Extra + s_ExtraSize, m_Extra);
    }

  ////////////////
  // operations //
  ////////////////

    //! cut out a part of the density map
    //! @param POSCOL   column or x position
    //! @param POSROW   row    pr y position
    //! @param POSLAYER layer  or z position
    //! @param NCOL     number of columns
    //! @param NROW     number of rows
    //! @param NLAYER   number of layers
    //! @return sub density map from given position of given size
    Map Map::SubMap
    (
      const size_t POSCOL, const size_t POSROW, const size_t POSLAYER,
      const size_t NCOL, const size_t NROW, const size_t NLAYER
    ) const
    {
      //make sure that coordinates does not extent the Denistymap
      BCL_Assert
      (
           POSCOL   < m_Data.GetNumberCols()   && POSCOL   + NCOL   <= m_Data.GetNumberCols()
        && POSROW   < m_Data.GetNumberRows()   && POSROW   + NROW   <= m_Data.GetNumberRows()
        && POSLAYER < m_Data.NumberLayers() && POSLAYER + NLAYER <= m_Data.NumberLayers(),
        "Sub coordinates extent Densitymap"
      );

      // extension output map
      const storage::VectorND< 3, size_t> ext_output( NCOL, NROW, NLAYER);

      // index changes
      linal::VectorND< int, 3> new_index( 0, 0, 0);
      new_index( 0) = m_Index( 0) + POSCOL;
      new_index( 1) = m_Index( 1) + POSROW;
      new_index( 2) = m_Index( 2) + POSLAYER;

      // intervals
      linal::VectorND< int, 3> intervals( 0, 0, 0);
      for( size_t i( 0); i < 3; ++i)
      {
        if( new_index( i) <= -int( ext_output( i)))
        {
          intervals( i) = math::Absolute( new_index( i));
        }
        else if( new_index( i) >= 0)
        {
          intervals( i) = new_index( i) + ext_output( i) - 1;
        }
        else
        {
          intervals( i) = ext_output( i) - 1;
        }
      }

      // length changes
      linal::Vector3D new_length( 0.0);
      new_length( 0) = intervals( 0) * m_CellWidth( 0);
      new_length( 1) = intervals( 1) * m_CellWidth( 1);
      new_length( 2) = intervals( 2) * m_CellWidth( 2);

      // create submap and return
      return Map
      (
        m_Data.SubTensor( NLAYER, NROW, NCOL, POSLAYER, POSROW, POSCOL),
        new_index,
        intervals,
        new_length,
        m_CellWidth,
        m_Angle,
        m_Axis,
        m_Origin
      );
    }

    //! get the common sub tensor from this and argument density map
    //! @param DENSITY_MAP the second density map to find common subtensor
    //! @return two tensors, first tensor of this map, second tensot of argument
    storage::VectorND< 2, math::Tensor< double> > Map::CommonSubTensor( const Map &DENSITY_MAP) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        math::EqualWithinTolerance( m_CellWidth, DENSITY_MAP.m_CellWidth),
        "Calculating common sub tensor requires density maps with equal CellWidth\nthis: " +
        util::Format()( m_CellWidth) + "\nargument: " + util::Format()( DENSITY_MAP.m_CellWidth)
      );

      // make sure that the origins are the same
      BCL_Assert
      (
        math::EqualWithinTolerance( m_Origin, DENSITY_MAP.m_Origin),
        "Calculating common sub tensor requires density maps with equal Origin\nthis: " +
        util::Format()( m_Origin) + "\nargument: " + util::Format()( DENSITY_MAP.m_Origin)
      );

      // relative index of argument density relative to this
      linal::Vector< int> this_start( 3, int( 0));
      linal::Vector< int> arg_start(  3, int( 0));
      linal::Vector< int> this_end(   3, int( 0));
      linal::Vector< int> arg_end(    3, int( 0));
      // find common start and end
      linal::Vector< int> common_start( 3, int( 0));
      linal::Vector< int> common_end( 3, int( 0));
      linal::Vector< int> extent( 3, int( 0));

      // iterate over dimensions
      for( size_t i( 0); i < 3; ++i)
      {
        this_start( i) =             m_Index( i);
        arg_start(  i) = DENSITY_MAP.m_Index( i);
        const int this_end      = this_start( i) +             GetDimensions()( i);
        const int arg_end       = arg_start(  i) + DENSITY_MAP.GetDimensions()( i);
        const int common_start  = std::max( this_start( i), arg_start( i));
        const int common_end    = std::min( this_end      , arg_end);
        extent( i) = common_end - common_start - 1;

        // no common sub density
        if( extent( i) <= 0)
        {
          return storage::VectorND< 2, math::Tensor< double> >();
        }

        this_start( i) = common_start -             m_Index( i);
        arg_start(  i) = common_start - DENSITY_MAP.m_Index( i);
      }

      // create sub tensors and return
      return storage::VectorND< 2, math::Tensor< double> >
      (
        m_Data.SubTensor( extent( 2), extent( 1), extent( 0), this_start( 2), this_start( 1), this_start( 0)),
        DENSITY_MAP.m_Data.SubTensor( extent( 2), extent( 1), extent( 0), arg_start( 2), arg_start( 1), arg_start( 0))
      );
    }

    //! @brief calculate the min max and mean from the intensities stored in the data
    void Map::CalculateMinMaxMeanRmsd()
    {
      m_Minimum = math::Statistics::MinimumValue( m_Data.Begin(), m_Data.End());
      m_Maximum = math::Statistics::MaximumValue( m_Data.Begin(), m_Data.End());
      m_Mean    = math::Statistics::Mean( m_Data.Begin(), m_Data.End());
      m_Rmsd    = math::Statistics::StandardDeviation( m_Data.Begin(), m_Data.End());
    }

    //! compute histogram over all densities in BIN bins
    math::Histogram Map::Histogram( const size_t NR_BINS) const
    {
      math::Histogram histogram( m_Minimum, ( m_Maximum - m_Minimum) / NR_BINS, NR_BINS);
      for( const double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        histogram.PushBack( *ptr);
      }

      return histogram;
    }

    //! compute a point cloud of NUMBER_OF_POINTS with an at least distance of FEATURE_DISTANCE in Angstrom using
    //! the highest derivatives between voxels in the map and the own intensity as the user can weight
    coord::PointCloud
    Map::CalculatePointCloud
    (
      const size_t NUMBER_OF_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT
    ) const
    {
      // local variables
      // function, that evaluates if a voxel overlaps by the given radius with another voxel
      const VoxelOverlap does_voxel_overlap_within_radius
      (
        FEATURE_DISTANCE,
        m_CellWidth
      );

      // set that keeps track of voxels with highest intensities and gradients
      std::multiset< Voxel, VoxelDensityGreater> considered_voxels;

      // set that stores the voxels that do overlap, but might be not considered too early
      std::multiset< Voxel, VoxelDensityGreater> overlapping_voxels;

      // store min intensity
      const double min_intensity( GetMean());

      // min intensity of a voxels has to be at least 50% of the highest intensity
      BCL_MessageVrb( "PointCloud( NUMBER_OF_POINTS = " + util::Format()( NUMBER_OF_POINTS) +
        ", FEATURE_DISTANCE = " + util::Format()( FEATURE_DISTANCE) + ")\nstart of building Pointcloud");

      // loop over all voxels and memorize best NUMBER_OF_POINTS
      for( size_t i( 1); i < m_Data.NumberLayers() - 1; ++i)
      {
        for( size_t j( 1); j < m_Data.GetNumberRows() - 1; ++j)
        {
          for( size_t k( 1); k < m_Data.GetNumberCols() - 1; ++k)
          {
            const double current_intensity( m_Data( i, j, k));

            if( current_intensity < min_intensity)
            {
              continue;
            }

            // make new voxel
            const double derivative_intensity
            (
              RATIO_INTENSITY_GRADIENT * current_intensity +
              double( 1) / double( 24) *
              (
                // direct nighbors
                math::Absolute( current_intensity - m_Data( i - 1, j    , k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k    )) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j    , k + 1))
              ) +
              double( 1) / ( double( 24) * math::Sqrt( double( 2))) *
              (
                // neighbors on edges
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k    )) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k    )) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i    , j + 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j    , k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j    , k + 1))
              ) +
              double( 1) / ( double( 24) * math::Sqrt( double( 3))) *
              (
                // neighbor on corners
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i - 1, j + 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j - 1, k + 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k - 1)) +
                math::Absolute( current_intensity - m_Data( i + 1, j + 1, k + 1))
              )
            );
            Voxel current_voxel = { i, j, k, derivative_intensity};

            std::vector< std::multiset< Voxel, VoxelDensityGreater>::iterator> current_overlaps;
            // find all overlapping voxels in the set of considered voxels
            std::multiset< Voxel, VoxelDensityGreater>::iterator
              voxel_itr_overlap( considered_voxels.begin()), voxel_itr_end( considered_voxels.end());
            bool all_smaller_intensity( true);
            while( true)
            {
              voxel_itr_overlap =
                std::find_if
                (
                  voxel_itr_overlap,
                  voxel_itr_end,
                  std::bind2nd( does_voxel_overlap_within_radius, current_voxel)
                );
              if( voxel_itr_overlap == voxel_itr_end)
              {
                break;
              }
              else
              {
                current_overlaps.push_back( voxel_itr_overlap);
                all_smaller_intensity &= voxel_itr_overlap->m_Intensity < current_voxel.m_Intensity;
              }
              ++voxel_itr_overlap;
            }

            // insert if no overlap was found
            if( current_overlaps.empty())
            {
              considered_voxels.insert( current_voxel);
            }
            else
            {
              // replace with the higher intensity
              if( all_smaller_intensity)
              {
                considered_voxels.insert( current_voxel);
                for
                (
                  std::vector< std::multiset< Voxel, VoxelDensityGreater>::iterator>::const_iterator
                    cur_overlap_itr( current_overlaps.begin()), cur_overlap_itr_end( current_overlaps.end());
                  cur_overlap_itr != cur_overlap_itr_end;
                  ++cur_overlap_itr
                )
                {
                  overlapping_voxels.insert( **cur_overlap_itr);
                  considered_voxels.erase( *cur_overlap_itr);
                }
              }
              // insert it in the overlapping voxels
              else
              {
                overlapping_voxels.insert( current_voxel);
              }
            }
            // check if considered voxels exceeds the considered size
            if( considered_voxels.size() > NUMBER_OF_POINTS)
            {
              std::multiset< Voxel, VoxelDensityGreater>::iterator itr( considered_voxels.end());
              --itr;
              overlapping_voxels.insert( *itr);
              considered_voxels.erase( itr);
            }
          }
          // remove overlapping voxels that will not be used any more
          if( !overlapping_voxels.empty() && !considered_voxels.empty())
          {
            // find first voxel in overlapping voxels that has higher intensity than the worst in the considered set
            std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_end( overlapping_voxels.end()), voxel_itr_overlap
            (
              std::find_if( overlapping_voxels.begin(), voxel_itr_end, std::bind1st( VoxelDensityGreater(), *( --considered_voxels.end())))
            );
            // erase the range with intensity that is too low
            overlapping_voxels.erase( voxel_itr_overlap, voxel_itr_end);
          }
        }
      }

      // if an overlapping Voxel has a higher intensity than one of the considered Voxels and non of the Voxels that
      // have a higher intensity than this, it has to be inserted and all overlapping Voxels with lower intensity have
      // to be erased

      // iterator on the overlapping voxels
      std::multiset< Voxel, VoxelDensityGreater>::iterator
        over_voxel_itr( overlapping_voxels.begin()), over_voxel_itr_end( overlapping_voxels.end());

      // iterator on the considered
      std::multiset< Voxel, VoxelDensityGreater>::iterator
        cons_voxel_itr_low( considered_voxels.begin());

      // iterate as long as there are overlapping voxels
      while( over_voxel_itr != over_voxel_itr_end)
      {
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_voxel_itr_end( considered_voxels.end());
        // find the first voxel with lower intensity in the considered voxels
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_voxel_itr
        (
          std::find_if( cons_voxel_itr_low, cons_voxel_itr_end, std::bind1st( VoxelDensityGreater(), *over_voxel_itr))
        );

        // no voxel with lower intensity
        if( cons_voxel_itr == cons_voxel_itr_end)
        {
          break;
        }

        // search in the range that has a higher intensity for an overlapping voxel
        std::multiset< Voxel, VoxelDensityGreater>::iterator cons_higher_overlapping_voxel_itr
        (
          std::find_if( considered_voxels.begin(), cons_voxel_itr, std::bind2nd( does_voxel_overlap_within_radius, *over_voxel_itr))
        );

        // if there is no overlapping iterator in that range
        if( cons_higher_overlapping_voxel_itr == cons_voxel_itr)
        {
          // insert that voxel
          considered_voxels.insert( *over_voxel_itr);
        }
        // proceed to next
        ++over_voxel_itr;
      }

      // delete all overlapping voxels with lower intensities
      std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr( considered_voxels.begin()), voxel_itr_end( considered_voxels.end());
      while( voxel_itr != voxel_itr_end)
      {
        std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_next( voxel_itr);
        ++voxel_itr_next;
        // at end
        if( voxel_itr_next == voxel_itr_end)
        {
          break;
        }

        // delete all overlapping points
        while( true)
        {
          std::multiset< Voxel, VoxelDensityGreater>::iterator voxel_itr_overlap
          (
            std::find_if( voxel_itr_next, voxel_itr_end, std::bind2nd( does_voxel_overlap_within_radius, *voxel_itr))
          );

          // remove overlapping
          if( voxel_itr_overlap != voxel_itr_end)
          {
            // progress after the overlapping iterator
            voxel_itr_next = voxel_itr_overlap;
            ++voxel_itr_next;
            considered_voxels.erase( voxel_itr_overlap);
          }
          // no remaining overlapping voxel
          else
          {
            break;
          }
        }
        ++voxel_itr;
      }

      BCL_MessageVrb
      (
        "NUMBER_OF_POINTS = " + util::Format()( NUMBER_OF_POINTS) + ", FEATURE_DISTANCE = " + util::Format()( FEATURE_DISTANCE) +
        "\nPoints found at all: " + util::Format()( considered_voxels.size()) + " will be reduced to target number of points"
      );

      // build PointCloud
      coord::PointCloud pointcloud;
      pointcloud.GetData().AllocateMemory( NUMBER_OF_POINTS);

      size_t count( 0);
      for
      (
        std::multiset< Voxel, VoxelDensityGreater>::const_iterator
          itr( considered_voxels.begin()), itr_end( considered_voxels.end());
        itr != itr_end && count < NUMBER_OF_POINTS;
        ++itr, ++count
      )
      {
        pointcloud.PushBack
        (
          linal::Vector3D
          (
            double( ( int( itr->m_Col)   + m_Index( 0)) - 0.5) * m_CellWidth.X() + m_Origin.X(),
            double( ( int( itr->m_Row)   + m_Index( 1)) - 0.5) * m_CellWidth.Y() + m_Origin.Y(),
            double( ( int( itr->m_Layer) + m_Index( 2)) - 0.5) * m_CellWidth.Z() + m_Origin.Z()
          )
        );
      }

      return pointcloud;
    }

    //! @brief add noise
    //! @param RNG an random number Generator object
    //! @param MEAN the mean of the gaussian distribution
    //! @param STANDARD_DEVIATION standard deviation of the gaussian distribution
    //! @return the cross correlation coefficient map (current map as simulated) vs new map with noise
    double Map::AddNoise( const random::DistributionInterface &RNG, const double MEAN, const double STANDARD_DEVIATION)
    {
      // copy of this map for cross correlation calculation after noise addition
      Map original_map( *this);

      // iterate over map and add random noise
      for( double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        *ptr += RNG.RandomGaussian( MEAN, STANDARD_DEVIATION);
      }

      // recalculate min, max, mean and rmsd
      CalculateMinMaxMeanRmsd();

      // end - calculate the correlation to the original map
      return CrossCorrelationCoefficient( original_map, original_map.m_Minimum);
    }

    //! @brief normalize
    //! transforms intensities, so that rmsd is 1
    void Map::Normalize()
    {
      // subtract mean from all intensities and devide by rmsd
      m_Data /= m_Rmsd;

      // recalculate min, max, mean and rmsd
      CalculateMinMaxMeanRmsd();
    }

    //! calculate the standard deviation to a smaller density map
    //! http://en.wikipedia.org/wiki/Cross-correlation#Normalized_cross-correlation
    double Map::Correlation( const Map &DENSITY_MAP) const
    {
      return CrossCorrelationCoefficient( DENSITY_MAP, double( 0));
    }

    //! calculate the cross correlation factor to a simulated density map, only using voxels that are above the given
    //! CONTOUR_LEVEL in the simulated density map
    //! @param SIMULATED_DENSITY_MAP
    //! @param CONTOUR_LEVEL
    //! @return cross correlation coefficient
    double Map::CrossCorrelationCoefficient
    (
      const Map &SIMULATED_DENSITY_MAP,
      const double CONTOUR_LEVEL
    ) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        math::EqualWithinTolerance( m_CellWidth, SIMULATED_DENSITY_MAP.m_CellWidth),
        "Calculating cross correlation factor requires density maps with equal CellWidth\nthis: " +
        util::Format()( m_CellWidth) + "\nargument: " + util::Format()( SIMULATED_DENSITY_MAP.m_CellWidth)
      );

      //store the index for the lower right and upper left corner of the argument density map
      linal::Vector< int> lower_right_corner_argument( 3, int( 0)), upper_left_corner_argument( 3, int( 0));

      for( size_t i( 0); i < 3; ++i)
      {
        lower_right_corner_argument( i) = int( ( SIMULATED_DENSITY_MAP.m_Origin( i) - m_Origin( i)) / m_CellWidth( i) + double( SIMULATED_DENSITY_MAP.m_Index( i) - m_Index( i)));
        upper_left_corner_argument( i) = lower_right_corner_argument( i) + SIMULATED_DENSITY_MAP.GetDimensions()( i) - 1;
      }

      // make sure that at all indices of any corner of the argument density map is within the Map
      if(
          !(
                upper_left_corner_argument(  0) >= 0 && lower_right_corner_argument( 0) <= int( GetDimensions()( 0))
             && upper_left_corner_argument(  1) >= 0 && lower_right_corner_argument( 1) <= int( GetDimensions()( 1))
             && upper_left_corner_argument(  2) >= 0 && lower_right_corner_argument( 2) <= int( GetDimensions()( 2))
            )
        )
      {
        return util::GetUndefined< double>();
      }

      // number of voxels that are above the CONTOUR_LEVEL in the experimental (this) density map
      size_t count_voxel( 0);

      // mean and standard deviation
      math::RunningAverageSD< double> mean_sd_this;
      math::RunningAverageSD< double> mean_sd_sim;

      // define the indices to iterate over
      const size_t i1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 0))));
      const size_t i2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
      const size_t i1_max( GetDimensions()( 0));
      const size_t i2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 0));

      const size_t j1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 1))));
      const size_t j2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
      const size_t j1_max( GetDimensions()( 1));
      const size_t j2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 1));

      const size_t k1_min( size_t( std::max( int( 0),  lower_right_corner_argument( 2))));
      const size_t k2_min( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
      const size_t k1_max( GetDimensions()( 2));
      const size_t k2_max( SIMULATED_DENSITY_MAP.GetDimensions()( 2));

      // calculate correlation by iterating over all density values
      for( size_t i1( i1_min), i2( i2_min); i1 < i1_max && i2 < i2_max; ++i1, ++i2)
      {
        for( size_t j1( j1_min), j2( j2_min); j1 < j1_max && j2 < j2_max; ++j1, ++j2)
        {
          for( size_t k1( k1_min), k2( k2_min); k1 < k1_max && k2 < k2_max; ++k1, ++k2)
          {
            const double this_int( operator()( i1, j1, k1));

            // skip exp intensities below contour level
            if( this_int < CONTOUR_LEVEL)
            {
              continue;
            }

            const double sim_int( SIMULATED_DENSITY_MAP( i2, j2, k2));
            ++count_voxel;
            mean_sd_this += this_int;
            mean_sd_sim += sim_int;
          }
        }
      }

      const double mean_this( mean_sd_this.GetAverage());
      const double mean_sim( mean_sd_sim.GetAverage());

      double correlation( 0);

      // calculate correlation by iterating over all density values
      for( size_t i1( i1_min), i2( i2_min); i1 < i1_max && i2 < i2_max; ++i1, ++i2)
      {
        for( size_t j1( j1_min), j2( j2_min); j1 < j1_max && j2 < j2_max; ++j1, ++j2)
        {
          for( size_t k1( k1_min), k2( k2_min); k1 < k1_max && k2 < k2_max; ++k1, ++k2)
          {
            const double this_int( operator()( i1, j1, k1));

            // skip exp intensities below contour level
            if( this_int < CONTOUR_LEVEL)
            {
              continue;
            }

            const double sim_int( SIMULATED_DENSITY_MAP( i2, j2, k2));
            correlation += ( this_int - mean_this) * ( sim_int - mean_sim);
          }
        }
      }

      // normalize by number of voxel
      correlation = double( 1) / ( double( count_voxel * mean_sd_this.GetStandardDeviation() * mean_sd_sim.GetStandardDeviation())) * correlation;

      return correlation;
    }

    //! calculate the standard deviation to a smaller density map
    double Map::StandardDeviation( const Map &DENSITY_MAP) const
    {
      // make sure that the two density maps have the same CellWidth
      BCL_Assert
      (
        m_CellWidth == DENSITY_MAP.m_CellWidth,
        "Calculating correlation factor requires density maps with equal CellWidth"
      );

      //store the index for the lower right and upper left corner of the argument density map
      linal::Vector< int> lower_right_corner_argument( 3, int( 0)), upper_left_corner_argument( 3, int( 0));
      for( size_t i( 0); i < 3; ++i)
      {
        lower_right_corner_argument( i) = int( ( DENSITY_MAP.m_Origin( i) - m_Origin( i)) / m_CellWidth( i) + double( int( DENSITY_MAP.m_Index( i)) - int( m_Index( i))) + 0.5);
        upper_left_corner_argument( i) = lower_right_corner_argument( i) + DENSITY_MAP.GetDimensions()( i) - 1;
      }

      // make sure that at least one index of any corner of the argument density map is within the Map
      if(
          !(
                upper_left_corner_argument(  0) >= 0 && lower_right_corner_argument( 0) <= int( GetDimensions()( 0))
             && upper_left_corner_argument(  1) >= 0 && lower_right_corner_argument( 1) <= int( GetDimensions()( 1))
             && upper_left_corner_argument(  2) >= 0 && lower_right_corner_argument( 2) <= int( GetDimensions()( 2))
            )
        )
      {
        return util::GetUndefined< double>();
      }

      std::pair< double, double> sim_min_max_threshold_intensity( DENSITY_MAP.GetMinimum(), DENSITY_MAP.GetMaximum());

      // variables to store average value and standard deviation of experimanetal and model map
      // but only for occupied voxels in simulatet map
      std::pair< double, double> exp_avg_sd( 0, 0), sim_avg_sd( 0, 0);
      size_t count_occupied_voxels( 0), count_voxels( 0);

      //calculate average density over all occupied voxels in simulated map
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
            {
              exp_avg_sd.first += operator()( i1, j1, k1);
              sim_avg_sd.first += DENSITY_MAP( i2, j2, k2);
              count_voxels++;
            }
      exp_avg_sd.first /= count_voxels;
      sim_avg_sd.first /= count_voxels;

      BCL_MessageVrb( "calculating standard deviation for " + util::Format()( count_voxels) + " overlapping voxels");

      //calculate standard deviation in occupied voxels
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
            {
              exp_avg_sd.second += math::Sqr( exp_avg_sd.first - operator()( i1, j1, k1));
              sim_avg_sd.second += math::Sqr( sim_avg_sd.first - DENSITY_MAP( i2, j2, k2));
            }
      //calculate average density in Map and Argument
      exp_avg_sd.second /= ( count_voxels - 1);
      sim_avg_sd.second /= ( count_voxels - 1);

      //calculate standard deviation in Map and Argument
      exp_avg_sd.first = math::Sqrt( exp_avg_sd.first);
      sim_avg_sd.first = math::Sqrt( sim_avg_sd.first);
      sim_min_max_threshold_intensity.first  = sim_avg_sd.first - 3 * sim_avg_sd.second;
      sim_min_max_threshold_intensity.second = sim_avg_sd.first + 3 * sim_avg_sd.second;

      double standard_deviation( 0);

      // calculate standard deviation
      for( size_t i1( size_t( std::max( int( 0), lower_right_corner_argument( 2)))),
                  i2( size_t( std::max( int( 0), -lower_right_corner_argument( 2))));
                  i1 < GetDimensions()( 2) && i2 < DENSITY_MAP.GetDimensions()( 2); ++i1, ++i2)
      {
        for( size_t j1( size_t( std::max( int( 0), lower_right_corner_argument( 1)))),
                    j2( size_t( std::max( int( 0), -lower_right_corner_argument( 1))));
                    j1 < GetDimensions()( 1) && j2 < DENSITY_MAP.GetDimensions()( 1); ++j1, ++j2)
        {
          for( size_t k1( size_t( std::max( int( 0), lower_right_corner_argument( 0)))),
                      k2( size_t( std::max( int( 0), -lower_right_corner_argument( 0))));
                      k1 < GetDimensions()( 0) && k2 < DENSITY_MAP.GetDimensions()( 0); ++k1, ++k2)
          {
            // shift gaussian distributions to same avarage value and scale to same standard deviation and substract values in voxels but only if target density map voxel is occupied
//            if( DENSITY_MAP( i2, j2, k2) > sim_min_max_threshold_intensity.first && DENSITY_MAP( i2, j2, k2) < sim_min_max_threshold_intensity.second)
            {
              standard_deviation += math::Sqr( ( ( DENSITY_MAP( i2, j2, k2) - sim_avg_sd.first) * exp_avg_sd.second) / sim_avg_sd.second + exp_avg_sd.first - operator()( i1, j1, k1));
              count_occupied_voxels++;
            }
          }
        }
      }

      standard_deviation = math::Sqrt( standard_deviation / count_occupied_voxels);

      return standard_deviation;
    }

    //! @brief determine edges using a Sobel operator
    //! @sa http://en.wikipedia.org/wiki/Edge_detection
    Map &Map::EdgeDetectionSobel( const double WEIGTH_FACE, const double WEIGHT_EDGE, const double WEIGHT_VERTEX)
    {
      math::Tensor< double> edges( m_Data);
      edges.SetZero();

      const double s_mask_x[ 27] =
      {
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         WEIGHT_EDGE,    WEIGTH_FACE,  WEIGHT_EDGE,
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,

         0,  0,  0,
         0,  0,  0,
         0,  0,  0,

        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX,
        -WEIGHT_EDGE,   -WEIGTH_FACE, -WEIGHT_EDGE,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX
      };

      const double s_mask_y[ 27] =
      {
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,
         WEIGHT_EDGE,    0, -WEIGHT_EDGE,
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,

         WEIGHT_EDGE,  0, -WEIGHT_EDGE,
         WEIGTH_FACE,  0, -WEIGTH_FACE,
         WEIGHT_EDGE,  0, -WEIGHT_EDGE,

         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX,
         WEIGHT_EDGE,    0, -WEIGHT_EDGE,
         WEIGHT_VERTEX,  0, -WEIGHT_VERTEX
      };

      const double s_mask_z[ 27] =
      {
         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         0, 0, 0,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX,

         WEIGHT_EDGE,  WEIGTH_FACE,  WEIGHT_EDGE,
         0, 0, 0,
        -WEIGHT_EDGE, -WEIGTH_FACE, -WEIGHT_EDGE,

         WEIGHT_VERTEX,  WEIGHT_EDGE,  WEIGHT_VERTEX,
         0, 0, 0,
        -WEIGHT_VERTEX, -WEIGHT_EDGE, -WEIGHT_VERTEX
      };

      const math::Tensor< double> mask_x( 3, 3, 3, s_mask_x);
      const math::Tensor< double> mask_y( 3, 3, 3, s_mask_y);
      const math::Tensor< double> mask_z( 3, 3, 3, s_mask_z);

      // iterate over all voxel
      for( size_t x = 1; x < m_Data.NumberLayers() - 1; ++x)
      {
        for( size_t y = 1; y < m_Data.GetNumberRows() - 1; ++y)
        {
          for( size_t z = 1; z < m_Data.GetNumberCols() - 1; ++z)
          {
            double sum_x( 0);
            double sum_y( 0);
            double sum_z( 0);

            for( int i = -1; i < 2; ++i)
            {
              for( int j = -1; j < 2; ++j)
              {
                for( int k = -1; k < 2; ++k)
                {
                  // x gradient
                  sum_x += m_Data( x + i, y + j, z + k) * mask_x( i + 1, j + 1, k + 1);
                  // y gradient
                  sum_y += m_Data( x + i, y + j, z + k) * mask_y( i + 1, j + 1, k + 1);
                  // z gradient
                  sum_z += m_Data( x + i, y + j, z + k) * mask_z( i + 1, j + 1, k + 1);
                }
              }
            }

            edges( x, y, z) = math::Absolute( sum_x) + math::Absolute( sum_y) + math::Absolute( sum_z);
          }
        }
      }

      // update the data
      m_Data = edges;

      // end
      return *this;
    }

    //! @brief balance the intensities
    //! @param LENGTH of cube vertex in which intensities are balanced
    Map &Map::BalanceIntensities( const double LENGTH)
    {
      if( LENGTH <= double( 0))
      {
        return *this;
      }

      const int nr_voxels_x( int( ( LENGTH / 2) / m_CellWidth.X()));
      const int nr_voxels_y( int( ( LENGTH / 2) / m_CellWidth.Y()));
      const int nr_voxels_z( int( ( LENGTH / 2) / m_CellWidth.Z()));

      if( nr_voxels_x == 0 || nr_voxels_y == 0 || nr_voxels_z == 0)
      {
        return *this;
      }
      BCL_MessageStd( "balancing intensities");
      math::Tensor< double> balanced( m_Data);
      const double size( ( 2 * nr_voxels_x + 1) * ( 2 * nr_voxels_y + 1) * ( 2 * nr_voxels_z + 1));
      // iterate over all voxel
      for( size_t x = 0; x < m_Data.GetNumberCols(); ++x)
      {
        for( size_t y = 0; y < m_Data.GetNumberRows(); ++y)
        {
          for( size_t z = 0; z < m_Data.NumberLayers(); ++z)
          {
            double sum( 0);
            for( int i = -nr_voxels_x; i < nr_voxels_x; ++i)
            {
              for( int j = -nr_voxels_y; j < nr_voxels_y; ++j)
              {
                for( int k = -nr_voxels_z; k < nr_voxels_z; ++k)
                {
                  sum += m_Data
                    (
                      ( i + x + m_Data.GetNumberCols()) % m_Data.GetNumberCols(),
                      ( j + y + m_Data.GetNumberRows()) % m_Data.GetNumberRows(),
                      ( k + z + m_Data.NumberLayers()) % m_Data.NumberLayers()
                    );
                }
              }
            }
            balanced( x, y, z) *= size / sum;
          }
        }
      }

      m_Data = balanced;
      BCL_MessageStd( "balancing intensities finished");

      // end
      return *this;
    }

    //! @brief convert Map into Spline
    math::TricubicSpline Map::ConvertDensityToSpline() const
    {
      const math::SplineBorderType border[ 3] = { math::e_FirstDer, math::e_FirstDer, math::e_FirstDer};

      // these vectors are used to input the starting points and the
      // grid width delta of every dimension (x, y, z) into the spline

      //calculate temporary origin
      const linal::Vector3D origin
      (
        m_Origin + linal::Vector3D
        (
          ( m_Index( 0) + 0.5) * m_CellWidth.X(),
          ( m_Index( 1) + 0.5) * m_CellWidth.Y(),
          ( m_Index( 2) + 0.5) * m_CellWidth.Z()
        )
      );
      const double start[] =
      {
        origin.Z(), origin.Y(), origin.X()
      };

      const double delta[] =
      {
        m_CellWidth.Z(), m_CellWidth.Y(), m_CellWidth.X()
      };

      const bool lin_cont[ 3] = { true, true, true};

      //this vector controls the behavior of the spline at the beginning and
      //end of every dimension, only has impact for SplineBorderType FIRSTDER

      //every pair describes the value of the first order derivative at start and end
      const storage::Pair< double, double> first_be[ 3] =
      {
        storage::Pair< double, double>( 1, -1),
        storage::Pair< double, double>( 1, -1),
        storage::Pair< double, double>( 1, -1)
      };

      //convert parameter density map in TriCubicSpline
      math::TricubicSpline densitymap_as_spline;
      densitymap_as_spline.Train( border, start, delta, m_Data, lin_cont, first_be);

      return densitymap_as_spline;
    }

    //! @brief orthogonalize map if angles aren't 90,90,90
    Map Map::OrthogonalizeMap() const
    {
      // if angles are already 90,90,90
      if( m_Angle.Z() == 90 && m_Angle.Y() == 90 && m_Angle.X() == 90)
      {
        // return just the original density map
        return Map
        (
          m_Data, m_Index, m_Intervals, m_Length, m_CellWidth, m_Angle, m_Axis, m_Origin
        );
      }

      // if angles are not already 90,90,90 then orthogonalize the map
      // for this calculate the point that is origin + index of the original map in real space coordinates
      // (this is point in map that is closest to the origin)
      const linal::Vector3D begin_of_map_non_orthog_system
      (
        m_CellWidth.X() * m_Index( 0) + m_Origin.X(),
        m_CellWidth.Y() * m_Index( 1) + m_Origin.Y(),
        m_CellWidth.Z() * m_Index( 2) + m_Origin.Z()
      );

      // begin of the original map in orthogonal system is the same (per definition)
      const linal::Vector3D begin_of_map_orthog_system( begin_of_map_non_orthog_system);

      // calculate the end point of original map (point that is farthest away from begin point and origin)
      const linal::Vector3D end_of_map_non_orthog_system
      (
        m_CellWidth.X() * ( m_Index( 0) + GetDimensions()( 0)) + m_Origin.X(),
        m_CellWidth.Y() * ( m_Index( 1) + GetDimensions()( 1)) + m_Origin.Y(),
        m_CellWidth.Z() * ( m_Index( 2) + GetDimensions()( 2)) + m_Origin.Z()
      );

      // convert the end point coordinates of original map into orthogonal system
      const linal::Vector3D end_of_map_orthog_system
      (
        end_of_map_non_orthog_system.X() + cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2)
          + cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1),
        end_of_map_non_orthog_system.Y() + cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2)
          + ( sin( math::Angle::Radian( m_Angle.X())) - 1.0) * m_CellWidth.Y() * GetDimensions()( 1),
        end_of_map_non_orthog_system.Z() + ( sin( math::Angle::Radian( m_Angle.Y())) - 1.0) * m_CellWidth.Z() * GetDimensions()( 2)
          + ( sin( math::Angle::Radian( m_Angle.Z())) - 1.0) * m_CellWidth.Z() * GetDimensions()( 2)
      );

      // calculate the new cell width (make sure it is not larger than the original cell width)
      const linal::Vector3D new_cell_width
      (
        std::abs( sin( math::Angle::Radian( m_Angle.X()))) * m_CellWidth.X(),
        std::abs( sin( math::Angle::Radian( m_Angle.Y()))) * m_CellWidth.Y(),
        std::abs( sin( math::Angle::Radian( m_Angle.Z()))) * m_CellWidth.Z()
      );

      // for debug purposes give out the begin and end point of the original map
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_map_orthog_system.Z()) + "  1.00 38.08           N"
      );
      BCL_MessageDbg
      (
        "ATOM   9996  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_map_orthog_system.Z()) + "  1.00 38.08           N"
      );

      // determine the begin and end points of the new map in real space coordinates (these points will be different
      // if any angle is greater 90 degrees)
      double begin_of_new_map_X( begin_of_map_orthog_system.X());
      double begin_of_new_map_Y( begin_of_map_orthog_system.Y());
      double begin_of_new_map_Z( begin_of_map_orthog_system.Z());
      double end_of_new_map_X( end_of_map_orthog_system.X());
      double end_of_new_map_Y( end_of_map_orthog_system.Y());
      double end_of_new_map_Z( end_of_map_orthog_system.Z());

      // if one of the angles is greater than 90, than the new map will be larger in one dimension
      // if angle is greater 90
      if( m_Angle.X() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_X += cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_X -= cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * GetDimensions()( 1);
      }

      // if angle is greater 90
      if( m_Angle.Y() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_X += cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_X -= cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * GetDimensions()( 2);
      }

      // if angle is greater 90
      if( m_Angle.Z() > 90)
      {
        // new begin will be old begin minus length of density map in this direction * cos ( 180 - angle)
        begin_of_new_map_Y += cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2);

        // new end will be old end plus length of density map in this direction * cos ( 180 - angle)
        end_of_new_map_Y -= cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * GetDimensions()( 2);
      }

      // for the new orthogonalized map make sure that it has no origin and all the shift is contained in the index
      // calculate the new index
      int new_index_X( ( begin_of_new_map_X) / new_cell_width.X());
      if( new_index_X > 0)
      {
        new_index_X = floor( new_index_X);
      }
      else if( new_index_X < 0)
      {
        new_index_X = ceil( new_index_X);
      }

      // calculate the new index
      int new_index_Y( ( begin_of_new_map_Y) / new_cell_width.Y());
      if( new_index_Y > 0)
      {
        new_index_Y = floor( new_index_Y);
      }
      else if( new_index_Y < 0)
      {
        new_index_Y = ceil( new_index_Y);
      }

      // calculate the new index
      int new_index_Z( ( begin_of_new_map_Z) / new_cell_width.Z());
      if( new_index_Z > 0)
      {
        new_index_Z = floor( new_index_Z);
      }
      else if( new_index_Z < 0)
      {
        new_index_Z = ceil( new_index_Z);
      }

      // calculate new extension of density map (corresponding to number of cols for the new map)
      size_t new_extension_X( ceil( ( end_of_new_map_X - begin_of_new_map_X) / new_cell_width.X()));

      // calculate new extension of density map (corresponding to number of rows for the new map)
      size_t new_extension_Y( ceil( ( end_of_new_map_Y - begin_of_new_map_Y) / new_cell_width.Y()));

      // calculate new extension of density map (corresponding to number of layers for the new map)
      size_t new_extension_Z( ceil( ( end_of_new_map_Z - begin_of_new_map_Z) / new_cell_width.Z()));

      // set begin of new map
      linal::Vector3D begin_of_new_map
      (
        begin_of_new_map_X,
        begin_of_new_map_Y,
        begin_of_new_map_Z
      );
      // set end of new map
      linal::Vector3D end_of_new_map
      (
        end_of_new_map_X,
        end_of_new_map_Y,
        end_of_new_map_Z
      );

      // for debug purposes give out the begin and end point of the new map
      BCL_MessageDbg
      (
        "ATOM   9997  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( begin_of_new_map.Z()) + "  1.00 38.08           N"
      );
      BCL_MessageDbg
      (
        "ATOM   9998  N   ASP B   1     " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.X()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.Y()) +
        " " + util::Format().W( 7).FFP( 3).R()( end_of_new_map.Z()) + "  1.00 38.08           N"
      );

      // set new index and extension
      linal::VectorND< int, 3> new_index( new_index_X, new_index_Y, new_index_Z);
      storage::VectorND< 3, size_t> new_extension( new_extension_X, new_extension_Y, new_extension_Z);

      // set the new intervals correctly
      linal::VectorND< int, 3> new_intervals;
      for( size_t i( 0); i < 3; ++i)
      {
        if( new_index( i) <= -int( new_extension( i)))
        {
          new_intervals( i) = math::Absolute( new_index( i));
        }
        else if( new_index( i) >= 0)
        {
          new_intervals( i) = new_index( i) + new_extension( i) - 1;
        }
        else
        {
          new_intervals( i) = new_extension( i) - 1;
        }
      }

      // set the new length
      const linal::Vector3D new_length
      (
        new_cell_width.X() * double( new_intervals( 0)),
        new_cell_width.Y() * double( new_intervals( 1)),
        new_cell_width.Z() * double( new_intervals( 2))
      );

      // set the new origin to 0, 0, 0
      linal::Vector3D new_origin( 0.0, 0.0, 0.0);

      // set the new angles to 90, 90, 90
      linal::Vector3D new_angles( 90, 90, 90);

      // convert the original map into a spline
      math::TricubicSpline old_map_as_spline( ConvertDensityToSpline());

      // create tensor of appropriate size for new data, filled with the minimum value of the original map
      math::Tensor< double> new_data
      (
        new_extension_Z, // 2
        new_extension_Y, // 1
        new_extension_X, // 0
        m_Minimum
      );

      // store dimensions of the original map in a storage::VectorND
      linal::VectorND< int, 3> old_dimensions
      (
        int( m_Data.GetNumberCols()),
        int( m_Data.GetNumberRows()),
        int( m_Data.NumberLayers())
      );

      // iterate over the x direction
      for( size_t itr_X( 0); itr_X < new_extension_X; ++itr_X)
      {
        // calculate real space coordinates of this particular point (for X coordinate)
        double coordinates_orthogonal_system_X( ( new_index_X + int( itr_X) + 0.5) * new_cell_width.X());

        // iterate over the y direction
        for( size_t itr_Y( 0); itr_Y < new_extension_Y; ++itr_Y)
        {
          // calculate real space coordinates of this particular point (for Y coordinate)
          double coordinates_orthogonal_system_Y( ( new_index_Y + int( itr_Y) + 0.5) * new_cell_width.Y());

          // iterate over the z direction
          for( size_t itr_Z( 0); itr_Z < new_extension_Z; ++itr_Z)
          {
            // calculate real space coordinates of this particular point (for Y coordinate)
            double coordinates_orthogonal_system_Z( ( new_index_Z + int( itr_Z) + 0.5) * new_cell_width.Z());

            // calculate real space coordinates of this particular point
            const linal::Vector3D coordinates_orthogonal_system
            (
              coordinates_orthogonal_system_X,
              coordinates_orthogonal_system_Y,
              coordinates_orthogonal_system_Z
            );

            // convert these real space coordinates into the non-orthogonal system
            const linal::Vector3D coordinates_non_orthogonal_system
            (
              coordinates_orthogonal_system.X() - cos( math::Angle::Radian( m_Angle.X())) * m_CellWidth.Y() * itr_Y
                - cos( math::Angle::Radian( m_Angle.Y())) * m_CellWidth.Z() * itr_Z,
              coordinates_orthogonal_system.Y() - cos( math::Angle::Radian( m_Angle.Z())) * m_CellWidth.Z() * itr_Z
                - ( sin( math::Angle::Radian( m_Angle.X())) - 1.0) * m_CellWidth.Y() * itr_Y,
              coordinates_orthogonal_system.Z() - ( sin( math::Angle::Radian( m_Angle.Y())) - 1.0) * m_CellWidth.Z() * itr_Z
                - ( sin( math::Angle::Radian( m_Angle.Z())) - 1.0) * m_CellWidth.Z() * itr_Z
            );

            // calculate the voxel in the old map that this real space coordinate corresponds to
            const linal::VectorND< int, 3> old_map_lookup
            (
              (( coordinates_non_orthogonal_system.X() - m_Origin.X()) / m_CellWidth.X()) - m_Index( 0),
              (( coordinates_non_orthogonal_system.Y() - m_Origin.Y()) / m_CellWidth.Y()) - m_Index( 1),
              (( coordinates_non_orthogonal_system.Z() - m_Origin.Z()) / m_CellWidth.Z()) - m_Index( 2)
            );

            // if you are looking up a coordinate that corresponds to a defined voxel in the original map, look it up
            // if this coordinate is outside the original map, don't do anything, new_data already has the minimum
            // value stored by default
            if
            (
              old_map_lookup( 0) >= 0 && old_map_lookup( 0) < old_dimensions( 0) &&
              old_map_lookup( 1) >= 0 && old_map_lookup( 1) < old_dimensions( 1) &&
              old_map_lookup( 2) >= 0 && old_map_lookup( 2) < old_dimensions( 2)
            )
            {
              // lookup intensity in spline and insert into this density map
              new_data( itr_Z, itr_Y, itr_X) =
                old_map_as_spline.F
                (
                  coordinates_non_orthogonal_system.Z(),
                  coordinates_non_orthogonal_system.Y(),
                  coordinates_non_orthogonal_system.X()
                );

              // this is second possible implementation that is not based on the spline
              // it just looks up the intensity at this particular voxel, no averaging for off-center-voxel value is done
//              new_data( itr_Z, itr_Y, itr_X) = m_Data( old_map_lookup( 2), old_map_lookup( 1), old_map_lookup( 0));
            }
          }
        }
      }

      // return the orthogonalized density map
      return Map
      (
        new_data, new_index, new_intervals, new_length, new_cell_width, new_angles, m_Axis, new_origin
      );
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Map::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "hashmap class that reads in mrc files");
      serializer.AddInitializer
        (
         "index",
         "index of map",
         io::Serialization::GetAgent( &m_Index)
         );
      serializer.AddInitializer
        (
         "intervals",
         "intervals along each axis",
         io::Serialization::GetAgent( &m_Intervals)
         );
      serializer.AddInitializer
        (
         "length",
         "size of unit cell in angstrom",
         io::Serialization::GetAgent( &m_Length)
         );
      serializer.AddInitializer
        (
         "cell width",
         "size of each voxel in angstrom",
         io::Serialization::GetAgent( &m_CellWidth)
         );
      serializer.AddInitializer
        (
         "angle",
         "voxel angle in degrees",
         io::Serialization::GetAgent( &m_Angle)
         );
      serializer.AddInitializer
        (
         "axis",
         "column, row, and section axes",
         io::Serialization::GetAgent( &m_Axis)
         );
      serializer.AddInitializer
        (
         "minimum",
         "minimum intensity",
         io::Serialization::GetAgent( &m_Minimum)
         );
      serializer.AddInitializer
        (
         "maximum",
         "maximum intensity",
         io::Serialization::GetAgent( &m_Maximum)
         );
      serializer.AddInitializer
        (
         "mean",
         "mean intensity",
         io::Serialization::GetAgent( &m_Mean)
         );
      serializer.AddInitializer
        (
         "space group",
         "space group number (0 or 1)",
         io::Serialization::GetAgent( &m_SpaceGroup)
         );
      serializer.AddInitializer
        (
         "number bytes symmetry data",
         "number of bytes used for symmetry data (0 or 80)",
         io::Serialization::GetAgent( &m_NumberBytesSymmetryData)
         );
      // serializer.AddInitializer
      //   (
      //    "extra",
      //    "extra data, 0 by default",
      //    io::Serialization::GetAgent( &m_Extra)
      //    );
      serializer.AddInitializer
        (
         "rmsd",
         "rmsd of density to mean density value",
         io::Serialization::GetAgent( &m_Rmsd)
         );
      serializer.AddInitializer
        (
         "origin",
         "origin of map",
         io::Serialization::GetAgent( &m_Origin)
         );
      serializer.AddInitializer
        (
         "data",
         "density as a tensor that also contains the dimensions of the density map",
         io::Serialization::GetAgent( &m_Data)
         );
      serializer.AddInitializer
        (
         "machine stamp",
         "machine stamp in header",
         io::Serialization::GetAgent( &m_MachineStamp)
         );
      serializer.AddInitializer
        (
         "labels",
         "labels in header",
         io::Serialization::GetAgent( &m_Labels)
         );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param MAP_RHS right hand side density map to be copied
    //! @return Map the changed map
    Map &Map::operator=( const Map &MAP_RHS)
    {
      // same map
      if( &MAP_RHS == this)
      {
        return *this;
      }

      // assign members
      m_Index                   = MAP_RHS.m_Index;
      m_Intervals               = MAP_RHS.m_Intervals;
      m_Length                  = MAP_RHS.m_Length;
      m_CellWidth               = MAP_RHS.m_CellWidth;
      m_Angle                   = MAP_RHS.m_Angle;
      m_Axis                    = MAP_RHS.m_Axis;
      m_Minimum                 = MAP_RHS.m_Minimum;
      m_Maximum                 = MAP_RHS.m_Maximum;
      m_Mean                    = MAP_RHS.m_Mean;
      m_SpaceGroup              = MAP_RHS.m_SpaceGroup;
      m_NumberBytesSymmetryData = MAP_RHS.m_NumberBytesSymmetryData;
      m_Rmsd                    = MAP_RHS.m_Rmsd;
      m_Origin                  = MAP_RHS.m_Origin;
      m_Data                    = MAP_RHS.m_Data;
      m_MachineStamp            = MAP_RHS.m_MachineStamp;
      m_Labels                  = MAP_RHS.m_Labels;

      // copy extra
      std::copy( MAP_RHS.m_Extra, MAP_RHS.m_Extra + s_ExtraSize, m_Extra);

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! Write Header
    std::ostream &Map::WriteHeader( std::ostream &OSTREAM) const
    {
      OSTREAM << "Dimensions in x, y, z:\n" << GetDimensions() << '\n'; // Dimensions of the map in x, y ,z
      OSTREAM << "Index:\n"                 << m_Index      << '\n'; // index of map
      OSTREAM << "Intervals:\n"             << m_Intervals  << '\n'; // intervals along each axis
      OSTREAM << "Length x, y, z:\n"        << m_Length     << '\n'; // size of unitcell in Angstrom
      OSTREAM << "Angle  a, b, c:\n"        << m_Angle      << '\n'; // voxel angle in degrees
      OSTREAM << "Axis( fast - slow):\n"    << m_Axis       << '\n'; // col, row and section axis
      OSTREAM << "Minimum Intensity:\n"     << m_Minimum    << '\n'; // Minimal density value
      OSTREAM << "Maximum Intensity:\n"     << m_Maximum    << '\n'; // Maximal density value
      OSTREAM << "Mean Intensity:\n"        << m_Mean       << '\n'; // Mean density value
      OSTREAM << "SpaceGroup: "             << m_SpaceGroup << '\n';
      OSTREAM << "Number sym data bytes: "  << m_NumberBytesSymmetryData << '\n';
      OSTREAM << "extra data: "             << std::string( m_Extra, s_ExtraSize) << '\n';
      OSTREAM << "RMSD:\n"                  << m_Rmsd         << '\n'; // rms of density to mean density value
      OSTREAM << "Origin x, y, z:\n"        << m_Origin       << '\n'; // origin of map
      OSTREAM << "machine stamp:\n"         << m_MachineStamp << '\n'; // machine stamp
      OSTREAM << "labels\n"                 << m_Labels       << '\n'; // labels

      // end
      return OSTREAM;
    }

    //! read Map from mrc file
    //! http://www2.mrc-lmb.cam.ac.uk/image2000.html
    std::istream &Map::ReadMRC( std::istream &ISTREAM, const size_t EXTENDED_HEADER)
    {
      // get length of file:
      ISTREAM.seekg( 0, std::ios::end);
      const std::istream::pos_type file_length( ISTREAM.tellg());
      std::istream::pos_type file_current;
      ISTREAM.seekg( 0, std::ios::beg);

      // determine if it is necessary to swap bytes when reading chars
      const bool swap( CheckSwap( ISTREAM));

      linal::VectorND< int, 3> dimensions;
      //read dimensions; bit 0-11
      dimensions( 0) = ReadInt( ISTREAM, swap); // NX
      dimensions( 1) = ReadInt( ISTREAM, swap); // NY
      dimensions( 2) = ReadInt( ISTREAM, swap); // NZ

      linal::VectorND< size_t, 3> dimensions_copy( dimensions.Begin(), dimensions.End());

      int mode; // MODE; bit 12-15
      mode = ReadInt( ISTREAM, swap);  //read mode of mrc file ( MODE     data type :
                                       //0 image : signed 8-bit bytes range -128 to 127
                                       //1 image : 16-bit halfwords
                                       //2 image : 32-bit reals
                                       //3 transform : complex 16-bit integers
                                       //4 transform : complex 32-bit reals

      BCL_Assert
      (
        mode == 2,
        "MRC file not written in MODE 2 - BCL is not able to read such a format. Written in MODE " +
        util::Format()( mode) + "\ndimensions read so far:\n" + util::Format()( dimensions) +
        "\nand swap is determined to be: " + util::Format()( swap)
      );

      // read Index; bit 16-27
      m_Index( 0)  = ReadInt( ISTREAM, swap); // NXSTART
      m_Index( 1) = ReadInt( ISTREAM, swap); // NYSTART
      m_Index( 2)  = ReadInt( ISTREAM, swap); // NZSTART

      // read Intervals; bit 28-39
      m_Intervals( 0)  = ReadInt( ISTREAM, swap); // MX
      m_Intervals( 1) = ReadInt( ISTREAM, swap); // MY
      m_Intervals( 2)  = ReadInt( ISTREAM, swap); // MZ

      //read width; bit bit 40-51
      for( double *ptr( m_Length.Begin()), *ptr_end( m_Length.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // CELLA
      }

      //calculate m_CellWidth
      m_CellWidth.X() = m_Length.X() / m_Intervals( 0);
      m_CellWidth.Y() = m_Length.Y() / m_Intervals( 1);
      m_CellWidth.Z() = m_Length.Z() / m_Intervals( 2);

      //read Angle; bit 52 - 63
      for( double *ptr( m_Angle.Begin()), *ptr_end( m_Angle.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // CELLB
      }

      //read Axis; bit 64-75
      m_Axis( 0)  = ReadInt( ISTREAM, swap); // MAPC cols
      m_Axis( 1) = ReadInt( ISTREAM, swap); // MAPC rows
      m_Axis( 2)  = ReadInt( ISTREAM, swap); // MAPC sections

      //read minimal, maximal and mean density value of map / will be recalculated; bit 76-87
      m_Minimum = ReadFloat( ISTREAM, swap); // DMIN
      m_Maximum = ReadFloat( ISTREAM, swap); // DMAX
      m_Mean    = ReadFloat( ISTREAM, swap); // DMEAN

      // additional information
      m_SpaceGroup = ReadInt( ISTREAM, swap); // ISPG; bit 88-91
      m_NumberBytesSymmetryData = ReadInt( ISTREAM, swap); // NSYMBT; bit 92-95
      BCL_Assert( m_NumberBytesSymmetryData == 0 || m_NumberBytesSymmetryData == 80, "number of bytes for symmetry data should either be 0 or 80");
      ISTREAM.read( m_Extra, s_ExtraSize); // bit 96-195

      // read origin of map; bit 196-207
      for( double *ptr( m_Origin.Begin()), *ptr_end( m_Origin.End()); ptr != ptr_end; ++ptr)
      {
        *ptr = ReadFloat( ISTREAM, swap); // ORIGIN
      }

      // read char string MAP; bit 208-211
      char tmp[ 5];
      tmp[ 4] = '\0';
      const char *map_ident( "MAP ");
      for( int i( 0); i < 4; ++i)
      {
        *( tmp + i) = ISTREAM.get();
      }
      if( std::string( tmp).compare( map_ident) != 0)
      {
        BCL_MessageCrt( "given mrc map does not contain \"MAP\" at pos 53 in header but: " + std::string( tmp));
      }

      // read machine stamp; bit 212-215
      m_MachineStamp = ReadInt( ISTREAM, swap);

      //checks if m_Axis is a right handed coordinate system - and orders if necessary
      if( m_Axis( 0) != 1 || m_Axis( 1) != 2)
      {
        BCL_MessageCrt
        (
          "Axis have a different order than x, y, z - all information will be ordered x, y, z. No Guarantee that this works!"
        );

        // apparently intervals, length , cellwidth and angle don't have to be ordered along with the axes
        // the index definitely has to be ordered
        // origin probably has to be ordered, but example that we used contained origin (0, 0, 0), so we couldn't
        // check whether origin has to be ordered...
        dimensions  = OrderValues( m_Axis, dimensions);
        m_Index     = OrderValues( m_Axis, m_Index);
//        m_Intervals = OrderValues( m_Axis, m_Intervals); // already ordered
//        m_Length    = OrderValues( m_Axis, m_Length);    // already ordered
//        m_CellWidth = OrderValues( m_Axis, m_CellWidth); // already ordered
//        m_Angle     = OrderValues( m_Axis, m_Angle);     // already ordered
        m_Origin    = OrderValues( m_Axis, m_Origin);
      }

      m_Data = math::Tensor< double>( dimensions( 2), dimensions( 1), dimensions( 0));

      //read rmsd; bit 216-219
      m_Rmsd = ReadFloat( ISTREAM, swap);

      // read number labels
      int number_labels;
      number_labels = ReadInt( ISTREAM, swap); // NLABL; bit 220-223
      BCL_Assert( number_labels <= s_MaxNumberLabels, "there should not be more than 10 labels, but there are: " + util::Format()( number_labels));

      // remove all current labels
      m_Labels.Reset();
      // read all labels
      for( int i( 0); i < number_labels; ++i)
      {
        char label[ s_LabelLength];
        ISTREAM.read( label, s_LabelLength);
        m_Labels.PushBack( std::string( label, s_LabelLength));
      }

      // read density into the tensor
      ISTREAM.seekg( 4 * 256 + EXTENDED_HEADER); //jump to end of header
      file_current = ISTREAM.tellg();

      if( int( file_length - file_current) != int( 4 * m_Data.GetSize()))
      {
        const std::istream::pos_type additional_header( ( file_length - file_current) - int( 4 * m_Data.GetSize()));
        std::ostringstream msg;
        msg << "header does not have a size of 1024 bits!\n"
            << "size of extended header is: " << additional_header << '\n'
            << "proceeding at logical position according to header information for number of stored intensities = "
            << GetSize() << '\n';
        BCL_MessageVrb( msg.str());
        ISTREAM.seekg( file_current + additional_header);
      }

      //read all intensity values
      if( m_Axis( 0) != 1 || m_Axis( 1) != 2)
      {
        linal::VectorND< size_t, 3> index;
        for( index( 2) = 0; index( 2) < dimensions_copy( 2); ++index( 2))
        {
          for( index( 1) = 0; index( 1) < dimensions_copy( 1); ++index( 1))
          {
            for( index( 0) = 0; index( 0) < dimensions_copy( 0); ++index( 0))
            {
              operator()( OrderValues( m_Axis, index)) = ReadFloat( ISTREAM, swap);
            }
          }
        }

        m_Axis = OrderValues( m_Axis, m_Axis);
      }
      else
      {
        for( double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end && !( ISTREAM.eof()); ++ptr)
        {
          *ptr = ReadFloat( ISTREAM, swap);
        }
      }

      // calculate the min max mean and rmsd
      CalculateMinMaxMeanRmsd();

      // only use index if origin was given, but index is 0
      if( m_Index( 0) == 0 && m_Index( 1) == 0 && m_Index( 2) == 0)
      {
        m_Index( 0) = int( m_Origin( 0) / m_CellWidth( 0));
        m_Index( 1) = int( m_Origin( 1) / m_CellWidth( 1));
        m_Index( 2) = int( m_Origin( 2) / m_CellWidth( 2));
        m_Origin = 0.0;
      }
      else // check if index is given, that origin is not used
      {
        if( !math::EqualWithinTolerance( 0.0, m_Origin.SquareNorm()))
        {
          BCL_MessageCrt
          (
            "map read that contains an origin and an index:\n" + util::Format()( m_Origin) + "\n" +
            util::Format()( m_Index) + "\norigin will be set to zero, assuming that the origin = index * voxel size"
          );
          m_Origin = 0.0;
        }
      }

      // end
      return ISTREAM;
    }

    //! write Map to mrc file
    //! http://www2.mrc-lmb.cam.ac.uk/image2000.html
    std::ostream &Map::WriteMRC( std::ostream &OSTREAM) const
    {
      int tmpint;
      float tmpfloat;
      tmpint = m_Data.GetNumberCols();       //write length in x direction;
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);
      tmpint = m_Data.GetNumberRows();       //write length in y direction
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);
      tmpint = m_Data.NumberLayers();     //write length in z direction;
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);

      tmpint = 2;
      //write mode of mrc file MODE     data type :
      //0 image : signed 8-bit bytes range -128 to 127
      //1 image : 16-bit halfwords
      //2 image : 32-bit reals
      //3 transform : complex 16-bit integers
      //4 transform : complex 32-bit reals
      OSTREAM.write( reinterpret_cast< const char *>( &tmpint), 4);

      //write Index
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 0)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Index( 2)), 4);

      //write Intervals
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 0)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Intervals( 2)), 4);

      //write width
      for( const double *ptr( m_Length.Begin()), *ptr_end( m_Length.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      //write Angle
      for( const double *ptr( m_Angle.Begin()), *ptr_end( m_Angle.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      //write Axis
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 0)),  4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 1)), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_Axis( 2)),  4);

      //write minimal, maximal and mean density value of map
      tmpfloat = float( m_Minimum);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      tmpfloat = float( m_Maximum);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      tmpfloat = float( m_Mean);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);

      // additional information
      OSTREAM.write( reinterpret_cast< const char *>( &m_SpaceGroup), 4);
      OSTREAM.write( reinterpret_cast< const char *>( &m_NumberBytesSymmetryData), 4);
      OSTREAM.write( m_Extra, s_ExtraSize);

      //write origin of map
      for( const double *ptr( m_Origin.Begin()), *ptr_end( m_Origin.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      // insert MAP
      const char *map_ident( "MAP ");
      OSTREAM.write( map_ident, 4);

      // insert machine stamp
      OSTREAM.write( reinterpret_cast< const char *>( &m_MachineStamp), 4);

      // write rms
      tmpfloat = float( m_Rmsd);
      OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);

      // write number labels
      const int number_labels( m_Labels.GetSize());
      OSTREAM.write( reinterpret_cast< const char *>( &number_labels), 4);

      // write all labels
      for( int i( 0); i < number_labels; ++i)
      {
        BCL_Assert( int( m_Labels( i).size()) == s_LabelLength, "label dose not have the required length!");
        OSTREAM.write( m_Labels( i).c_str(), s_LabelLength);
      }

      // fill stream with empty labels
      OSTREAM.write( std::string( ( s_MaxNumberLabels - number_labels) * s_LabelLength, ' ').c_str(), ( s_MaxNumberLabels - number_labels) * s_LabelLength);

      // write density from the tensor
      for( const double *ptr( m_Data.Begin()), *ptr_end( m_Data.End()); ptr != ptr_end; ++ptr)
      {
        tmpfloat = float( *ptr);
        OSTREAM.write( reinterpret_cast< const char *>( &tmpfloat), 4);
      }

      // end
      return OSTREAM;
    }

    //! write Map to std::ostream using the given util::Format
    std::ostream &Map::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Index,     OSTREAM, INDENT) << '\n'; // index of map
      io::Serialize::Write( m_Intervals, OSTREAM, INDENT) << '\n'; // intervals along each axis
      io::Serialize::Write( m_Length,    OSTREAM, INDENT) << '\n'; // size of unitcell in Angstrom
      io::Serialize::Write( m_CellWidth, OSTREAM, INDENT) << '\n'; // size of each voxel in Angstroem
      io::Serialize::Write( m_Angle,     OSTREAM, INDENT) << '\n'; // voxel angle in degrees
      io::Serialize::Write( m_Axis,      OSTREAM, INDENT) << '\n'; // col, row and section axis
      io::Serialize::Write( m_Minimum,   OSTREAM, INDENT) << '\n'; // Minimal density value
      io::Serialize::Write( m_Maximum,   OSTREAM, INDENT) << '\n'; // Maximal density value
      io::Serialize::Write( m_Mean,      OSTREAM, INDENT) << '\n'; // Mean density value
      io::Serialize::Write( m_Rmsd,      OSTREAM, INDENT) << '\n'; // rmsd of density to mean density value
      io::Serialize::Write( m_Origin,    OSTREAM, INDENT) << '\n'; // origin of map

      // density as a tensor does also contain the dimensions of the density map
      io::Serialize::Write( m_Data,      OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read Map from io::IFStream
    std::istream &Map::Read( std::istream &ISTREAM)
    {
      //read members
      io::Serialize::Read( m_Index,     ISTREAM); // index of map
      io::Serialize::Read( m_Intervals, ISTREAM); // intervals along each axis
      io::Serialize::Read( m_Length,    ISTREAM); // size of unitcell in Angstrom
      io::Serialize::Read( m_CellWidth, ISTREAM); // size of each voxel in Angstroem
      io::Serialize::Read( m_Angle,     ISTREAM); // voxel angle in degrees
      io::Serialize::Read( m_Axis,      ISTREAM); // col, row and section axis
      io::Serialize::Read( m_Minimum,   ISTREAM); // Minimal density value
      io::Serialize::Read( m_Maximum,   ISTREAM); // Maximal density value
      io::Serialize::Read( m_Mean,      ISTREAM); // Mean density value
      io::Serialize::Read( m_Rmsd,      ISTREAM); // rmsd of density to mean density value
      io::Serialize::Read( m_Origin,    ISTREAM); // origin of map
      io::Serialize::Read( m_Data,      ISTREAM); // density as a tensor does also contain the dimensions of the density map

      // end
      return ISTREAM;
    }

    //!checks stream of mrcfile wether the bytes are swapped
    bool Map::CheckSwap( std::istream &ISTREAM)
    {
      unsigned char *cptr;
      int unswapped, swapped;

      //read first number in header
      ISTREAM.read( reinterpret_cast< char *>( &unswapped), 4);
      //set pointer on stream back to begin of stream
      ISTREAM.seekg( 0);
      swapped = unswapped;

      //swap byte 0, 3 and 1, 2
      cptr = ( unsigned char *)&swapped;
      std::swap( cptr[0], cptr[3]);
      std::swap( cptr[1], cptr[2]);

      BCL_MessageDbg
      (
        "unswapped char: " + util::Format()( unswapped) + " swapped char: " + util::Format()( swapped)
      );

      //if swapped is smaller unswapped then all bytes of the file have to be swapped
      return std::abs( ( long int)( swapped)) < std::abs( ( long int)( unswapped));
    }

    //!read int from istream and swap if necessary
    int Map::ReadInt( std::istream &ISTREAM, const bool SWAP)
    {
      unsigned char *cptr;
      int tmp;

      ISTREAM.read( reinterpret_cast< char *>( &tmp), 4);
      if( SWAP)
      {
        cptr = ( unsigned char *)&tmp;
        std::swap( cptr[0], cptr[3]);
        std::swap( cptr[1], cptr[2]);
      }
      return tmp;
    }

    //!read float from istream and swap if necessary
    float Map::ReadFloat( std::istream &ISTREAM, const bool SWAP)
    {
      unsigned char *cptr;
      float tmp;

      ISTREAM.read( reinterpret_cast< char *>( &tmp), 4);
      if( SWAP)
      {
        cptr = ( unsigned char *)&tmp;
        std::swap( cptr[0], cptr[3]);
        std::swap( cptr[1], cptr[2]);
      }
      return tmp;
    }

    //!order values according to given Axis order
    linal::Vector3D Map::OrderValues( const linal::VectorND< int, 3> &AXIS, const linal::Vector3D &VALUES)
    {
      storage::VectorND< 3, size_t> index_xyz;
      index_xyz( AXIS( 0) - 1) = 0;
      index_xyz( AXIS( 1) - 1) = 1;
      index_xyz( AXIS( 2) - 1) = 2;

      linal::Vector3D new_values( VALUES);
      new_values( 0) = VALUES( index_xyz( 0));
      new_values( 1) = VALUES( index_xyz( 1));
      new_values( 2) = VALUES( index_xyz( 2));

      return new_values;
    }

  } // namespace density
} // namespace bcl
