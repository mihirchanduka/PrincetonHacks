"""
Advanced Behavioral Simulator for Synthetic Mice

This module provides a comprehensive behavioral simulation system that models
realistic mouse behaviors based on actual behavioral data and phenotypes.
"""
import numpy as np
import pandas as pd
import random
from typing import Dict, List, Tuple
import math


class BehavioralSimulator:
    """
    A sophisticated behavioral simulator that models realistic mouse behaviors
    based on phenotype data and incorporating multiple behavioral parameters.
    """
    
    def __init__(self):
        # Define behavioral parameters based on real mouse data
        self.behavioral_parameters = {
            'locomotion': {
                'baseline_speed': 5.0,  # cm/s
                'anxiety_speed_factor': -2.0,  # anxious mice move slower
                'activity_speed_factor': 1.5,  # active mice move faster
                'max_speed': 20.0  # maximum speed (cm/s)
            },
            'exploration': {
                'center_aversion': 0.8,  # higher for more anxious mice
                'novelty_exploration': 0.7,  # willingness to explore new areas
                'thigmotaxis_strength': 0.9  # tendency to stay near walls
            },
            'social': {
                'social_preference': 0.6,  # baseline social behavior
                'aggression_factor': 0.8  # higher aggression affects social behavior
            }
        }
    
    def simulate_open_field_test(
        self, 
        phenotype: Dict[str, float], 
        duration: int = 300,  # 5 minutes
        arena_size: int = 100  # cm
    ) -> List[Dict[str, float]]:
        """
        Simulate an open field test with realistic behavioral patterns.
        
        Args:
            phenotype: Dictionary of phenotype scores
            duration: Duration in seconds
            arena_size: Arena size in cm
            
        Returns:
            List of behavioral data points
        """
        # Extract phenotype scores
        anxiety = phenotype.get('anxiety_score', phenotype.get('anxiety', 0.5))
        memory = phenotype.get('memory_score', phenotype.get('memory', 0.5))
        activity = phenotype.get('activity_score', phenotype.get('activity', 0.5))
        
        # Initialize position
        x, y = arena_size / 2, arena_size / 2  # Start in center
        center_zone_boundary = arena_size / 4  # Center zone is 1/4 of arena
        
        # Behavioral parameters based on phenotype
        speed_factor = (1 - anxiety) * self.behavioral_parameters['locomotion']['anxiety_speed_factor'] + \
                      activity * self.behavioral_parameters['locomotion']['activity_speed_factor']
        
        # Calculate behavioral tendencies
        thigmotaxis = self.behavioral_parameters['exploration']['thigmotaxis_strength'] * anxiety
        center_aversion = self.behavioral_parameters['exploration']['center_aversion'] * anxiety
        exploration_drive = self.behavioral_parameters['exploration']['novelty_exploration'] * (1 - anxiety)
        
        behavior_data = []
        velocity_x, velocity_y = 0.0, 0.0  # Initial velocity
        
        # Track visited zones to model habituation
        visited_zones = np.zeros((int(arena_size/10), int(arena_size/10)))
        habituation_map = np.ones((int(arena_size/10), int(arena_size/10)))
        
        for t in range(duration):
            # Determine current zone
            is_in_center = (center_zone_boundary < x < arena_size - center_zone_boundary and
                          center_zone_boundary < y < arena_size - center_zone_boundary)
            zone = "center" if is_in_center else "periphery"
            
            # Update habituation based on current location
            zone_x = min(int(x / 10), visited_zones.shape[0] - 1)
            zone_y = min(int(y / 10), visited_zones.shape[1] - 1)
            visited_zones[zone_x, zone_y] += 1
            # Habituation makes mice less likely to revisit the same area
            habituation_map[zone_x, zone_y] = max(0.3, 1.0 - visited_zones[zone_x, zone_y] * 0.1)
            
            # Calculate current speed based on factors
            base_speed = self.behavioral_parameters['locomotion']['baseline_speed']
            current_speed = max(0.5, base_speed + speed_factor) * habituation_map[zone_x, zone_y]
            
            # Movement behavior based on anxiety and exploration
            if is_in_center:
                # Mouse in center - affected by anxiety
                if random.random() < center_aversion:
                    # Move towards periphery (wall-hugging)
                    target_x, target_y = arena_size / 2, arena_size / 2  # Move toward center initially
                    if x < target_x:
                        target_x = arena_size * 0.9  # Move to right edge
                    else:
                        target_x = arena_size * 0.1  # Move to left edge
                    if y < target_y:
                        target_y = arena_size * 0.9  # Move to bottom edge
                    else:
                        target_y = arena_size * 0.1  # Move to top edge
                    
                    # Add some random movement
                    target_x += random.uniform(-5, 5)
                    target_y += random.uniform(-5, 5)
                    
                    # Calculate direction vector
                    dx = target_x - x
                    dy = target_y - y
                    dist = max(0.1, math.sqrt(dx*dx + dy*dy))
                    
                    # Normalize and scale by speed
                    velocity_x = (dx / dist) * current_speed
                    velocity_y = (dy / dist) * current_speed
                else:
                    # Exploratory movement in center
                    velocity_x += random.uniform(-2, 2) * exploration_drive
                    velocity_y += random.uniform(-2, 2) * exploration_drive
                    
                    # Limit speed
                    vel_mag = max(0.1, math.sqrt(velocity_x**2 + velocity_y**2))
                    if vel_mag > current_speed:
                        velocity_x = (velocity_x / vel_mag) * current_speed
                        velocity_y = (velocity_y / vel_mag) * current_speed
            else:
                # Mouse in periphery - thigmotaxis (wall-following) behavior
                if random.random() < thigmotaxis:
                    # Follow walls with some randomness
                    if x < center_zone_boundary:  # Left edge
                        velocity_x = max(0, random.uniform(-1, 2))
                        velocity_y = random.uniform(-2, 2)
                    elif x > arena_size - center_zone_boundary:  # Right edge
                        velocity_x = min(0, random.uniform(-2, 1))
                        velocity_y = random.uniform(-2, 2)
                    elif y < center_zone_boundary:  # Top edge
                        velocity_x = random.uniform(-2, 2)
                        velocity_y = max(0, random.uniform(-1, 2))
                    elif y > arena_size - center_zone_boundary:  # Bottom edge
                        velocity_x = random.uniform(-2, 2)
                        velocity_y = min(0, random.uniform(-2, 1))
                    else:  # Somewhere in the middle but considered periphery
                        # Move toward nearest wall
                        dist_left = x
                        dist_right = arena_size - x
                        dist_top = y
                        dist_bottom = arena_size - y
                        min_dist = min(dist_left, dist_right, dist_top, dist_bottom)
                        
                        if min_dist == dist_left:
                            velocity_x = random.uniform(-2, 0)
                            velocity_y = random.uniform(-1, 1)
                        elif min_dist == dist_right:
                            velocity_x = random.uniform(0, 2)
                            velocity_y = random.uniform(-1, 1)
                        elif min_dist == dist_top:
                            velocity_x = random.uniform(-1, 1)
                            velocity_y = random.uniform(-2, 0)
                        else:  # bottom
                            velocity_x = random.uniform(-1, 1)
                            velocity_y = random.uniform(0, 2)
                else:
                    # Move toward center occasionally
                    center_dx = (arena_size / 2) - x
                    center_dy = (arena_size / 2) - y
                    dist_to_center = max(0.1, math.sqrt(center_dx**2 + center_dy**2))
                    
                    velocity_x = (center_dx / dist_to_center) * current_speed
                    velocity_y = (center_dy / dist_to_center) * current_speed
            
            # Apply some random perturbation to make movement more natural
            velocity_x += random.uniform(-0.5, 0.5)
            velocity_y += random.uniform(-0.5, 0.5)
            
            # Limit maximum velocity
            vel_mag = math.sqrt(velocity_x**2 + velocity_y**2)
            if vel_mag > current_speed:
                velocity_x = (velocity_x / vel_mag) * current_speed
                velocity_y = (velocity_y / vel_mag) * current_speed
            
            # Update position
            new_x = x + velocity_x
            new_y = y + velocity_y
            
            # Keep within bounds with bouncing effect
            if new_x < 0 or new_x > arena_size:
                velocity_x *= -0.5  # Reverse and dampen
                new_x = max(0, min(arena_size, new_x))
            if new_y < 0 or new_y > arena_size:
                velocity_y *= -0.5  # Reverse and dampen
                new_y = max(0, min(arena_size, new_y))
            
            x, y = new_x, new_y
            
            # Record behavioral data point
            behavior_data.append({
                'x': round(x, 2),
                'y': round(y, 2),
                'time': t,
                'zone': zone,
                'velocity_x': round(velocity_x, 2),
                'velocity_y': round(velocity_y, 2),
                'speed': round(math.sqrt(velocity_x**2 + velocity_y**2), 2)
            })
        
        return behavior_data
    
    def simulate_elevated_plus_maze(
        self, 
        phenotype: Dict[str, float], 
        duration: int = 300
    ) -> List[Dict[str, any]]:
        """
        Simulate an elevated plus maze test - another common anxiety test.
        
        Args:
            phenotype: Dictionary of phenotype scores
            duration: Duration in seconds
            
        Returns:
            List of behavioral data points for elevated plus maze
        """
        anxiety = phenotype.get('anxiety_score', phenotype.get('anxiety', 0.5))
        exploration = 1 - anxiety  # Higher for less anxious mice
        
        # Maze structure: 4 arms (2 open, 2 closed)
        # For simplicity, we'll model movement between arm types
        arms = ['open_arm_1', 'open_arm_2', 'closed_arm_1', 'closed_arm_2']
        current_arm = 'closed_arm_1'  # Start in a closed arm (safer)
        
        behavior_data = []
        time_in_open = 0
        entries_to_open = 0
        
        for t in range(duration):
            # Decide movement based on anxiety
            # Anxious mice avoid open arms
            if 'open' in current_arm:
                time_in_open += 1
                # High anxiety makes it more likely to leave open arms quickly
                prob_stay_open = exploration
                if random.random() > prob_stay_open:
                    # Move to a closed arm
                    current_arm = random.choice(['closed_arm_1', 'closed_arm_2'])
            else:  # In closed arm
                # Lower anxiety makes it more likely to enter open arms
                prob_enter_open = exploration * 0.7  # Scaled down for realism
                if random.random() < prob_enter_open:
                    # Enter an open arm
                    current_arm = random.choice(['open_arm_1', 'open_arm_2'])
                    entries_to_open += 1
            
            behavior_data.append({
                'time': t,
                'arm_type': 'open' if 'open' in current_arm else 'closed',
                'arm': current_arm,
                'anxiety_indicator': 1 if 'open' in current_arm else 0  # Time in open arms
            })
        
        # Add summary metrics
        behavior_data.append({
            'summary': True,
            'time_in_open': time_in_open,
            'entries_to_open': entries_to_open,
            'open_time_ratio': time_in_open / duration if duration > 0 else 0
        })
        
        return behavior_data
    
    def simulate_novel_object_test(
        self, 
        phenotype: Dict[str, float], 
        duration: int = 300
    ) -> List[Dict[str, any]]:
        """
        Simulate a novel object test for assessing recognition memory.
        
        Args:
            phenotype: Dictionary of phenotype scores
            duration: Duration in seconds
            
        Returns:
            List of behavioral data points for novel object test
        """
        memory = phenotype.get('memory_score', phenotype.get('memory', 0.5))
        anxiety = phenotype.get('anxiety_score', phenotype.get('anxiety', 0.5))
        
        # Initialize positions
        arena_size = 100
        x, y = arena_size / 2, arena_size / 2  # Start in center
        center_zone_boundary = arena_size / 4
        
        # Place novel objects
        object_positions = [
            (arena_size * 0.25, arena_size * 0.25),  # Top-left
            (arena_size * 0.75, arena_size * 0.75)   # Bottom-right
        ]
        
        behavior_data = []
        object_interactions = [0, 0]  # Count interactions with each object
        
        for t in range(duration):
            is_in_center = (center_zone_boundary < x < arena_size - center_zone_boundary and
                          center_zone_boundary < y < arena_size - center_zone_boundary)
            zone = "center" if is_in_center else "periphery"
            
            # Check if close to an object (within 10 units)
            for i, obj_pos in enumerate(object_positions):
                dist_to_obj = math.sqrt((x - obj_pos[0])**2 + (y - obj_pos[1])**2)
                if dist_to_obj < 10:  # Within interaction zone
                    object_interactions[i] += 1
                    # Mouse tends to stay near object briefly
                    break
            else:
                # Not near any object, move according to anxiety/exploration
                if is_in_center and random.random() < (1 - anxiety):  # Less anxious explore more
                    # Move randomly in center to potentially encounter objects
                    x += random.uniform(-3, 3)
                    y += random.uniform(-3, 3)
                elif not is_in_center:
                    # High anxiety keeps mouse in periphery
                    if random.random() < anxiety:
                        # Move along periphery
                        if x < arena_size / 2:
                            x += random.uniform(-1, 2)
                        else:
                            x += random.uniform(-2, 1)
                        if y < arena_size / 2:
                            y += random.uniform(-1, 2)
                        else:
                            y += random.uniform(-2, 1)
                    else:
                        # Move toward center
                        x += (arena_size / 2 - x) * 0.1
                        y += (arena_size / 2 - y) * 0.1
            
            # Keep within bounds
            x = max(0, min(x, arena_size))
            y = max(0, min(y, arena_size))
            
            behavior_data.append({
                'x': round(x, 2),
                'y': round(y, 2),
                'time': t,
                'zone': zone,
                'object_0_interactions': object_interactions[0],
                'object_1_interactions': object_interactions[1],
                'memory_performance': memory  # Reference for analysis
            })
        
        return behavior_data
    
    def simulate_behavior(
        self, 
        phenotype: Dict[str, float], 
        duration: int = 300,
        test_type: str = "open_field"
    ) -> List[Dict[str, any]]:
        """
        General simulation function that can run different behavioral tests.
        
        Args:
            phenotype: Dictionary of phenotype scores
            duration: Duration in seconds
            test_type: Type of behavioral test to run
            
        Returns:
            List of behavioral data points
        """
        if test_type == "open_field":
            return self.simulate_open_field_test(phenotype, duration)
        elif test_type == "elevated_plus_maze":
            return self.simulate_elevated_plus_maze(phenotype, duration)
        elif test_type == "novel_object":
            return self.simulate_novel_object_test(phenotype, duration)
        else:
            # Default to open field
            return self.simulate_open_field_test(phenotype, duration)
    
    def calculate_behavioral_metrics(
        self, 
        behavior_data: List[Dict[str, any]], 
        test_type: str = "open_field"
    ) -> Dict[str, float]:
        """
        Calculate relevant behavioral metrics from the simulation data.
        
        Args:
            behavior_data: List of behavior data points
            test_type: Type of test to calculate metrics for
            
        Returns:
            Dictionary of behavioral metrics
        """
        if not behavior_data:
            return {}
        
        if test_type == "open_field":
            # Filter out any summary records
            actual_data = [d for d in behavior_data if 'summary' not in d]
            
            if not actual_data:
                return {}
            
            # Calculate metrics
            center_time = sum(1 for d in actual_data if d['zone'] == 'center')
            total_time = len(actual_data)
            center_time_ratio = center_time / total_time if total_time > 0 else 0
            
            # Calculate distance traveled
            total_distance = 0.0
            if len(actual_data) > 1:
                for i in range(1, len(actual_data)):
                    dx = actual_data[i]['x'] - actual_data[i-1]['x']
                    dy = actual_data[i]['y'] - actual_data[i-1]['y']
                    total_distance += math.sqrt(dx*dx + dy*dy)
            
            # Calculate average speed
            speeds = [d['speed'] for d in actual_data if 'speed' in d]
            avg_speed = sum(speeds) / len(speeds) if speeds else 0
            
            # Calculate zone transitions (indicates exploration)
            transitions = 0
            for i in range(1, len(actual_data)):
                if actual_data[i]['zone'] != actual_data[i-1]['zone']:
                    transitions += 1
            
            return {
                'center_time_ratio': center_time_ratio,
                'total_distance': round(total_distance, 2),
                'avg_speed': round(avg_speed, 2),
                'zone_transitions': transitions,
                'total_time': total_time
            }
        
        elif test_type == "elevated_plus_maze":
            # Find the summary record for this test
            summary_records = [d for d in behavior_data if d.get('summary', False)]
            if summary_records:
                summary = summary_records[0]
                return {
                    'time_in_open': summary['time_in_open'],
                    'entries_to_open': summary['entries_to_open'],
                    'open_time_ratio': summary['open_time_ratio']
                }
            return {}
        
        elif test_type == "novel_object":
            if behavior_data:
                final_record = behavior_data[-1]
                if 'object_0_interactions' in final_record and 'object_1_interactions' in final_record:
                    total_interactions = final_record['object_0_interactions'] + final_record['object_1_interactions']
                    return {
                        'object_0_interactions': final_record['object_0_interactions'],
                        'object_1_interactions': final_record['object_1_interactions'],
                        'total_interactions': total_interactions,
                        'memory_performance': final_record.get('memory_performance', 0.5)
                    }
            return {}
        
        return {}