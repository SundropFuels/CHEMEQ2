def parse_species(species_str):
    """Parses a molecule into into a dictionary of {element: # of atoms}"""

    parsed_elements = {}
    current_element = ""
    current_number_str = ""
    

    i = 0
    if not species_str[0].isalpha():
        raise BadCharacterError, "A molecule must start with an alphabetical character"

    while i < len(species_str):
        current_char = species_str[i]
        
        if current_char.isalpha():
            if i+1 == len(species_str) or species_str[i+1].isupper():
                #Allows for single character names, like CH4
                current_element = "".join([current_element, current_char])
                if current_element in parsed_elements.keys():
                    parsed_elements[current_element] += 1
                else:
                    parsed_elements[current_element] = 1
                current_element = ""
                
                
            else:
                
                current_element = "".join([current_element, current_char])
            i += 1
            continue

        elif current_char.isdigit():
            #we have gotten to the end of an element name
            
            while i < len(species_str) and species_str[i].isdigit():
                
                current_char = species_str[i]
                current_number_str = "".join([current_number_str, current_char])
                i += 1
            if current_number_str == '':
                raise BadCharacterError, "Each element must have a number associated with it"
            if current_element in parsed_elements.keys():
                parsed_elements[current_element] += int(current_number_str)
            else:
                parsed_elements[current_element] = int(current_number_str)
            
            current_element = ""
            current_number_str = ""
        else:
            raise BadCharacterError, "A molecule can only contain alphabetical and numerical characters"

    return parsed_elements


if __name__ == "__main__":
    f = parse_species("C2H4")

    print f

    g = parse_species("CH3OH")
    print g

    m = parse_species("Ar")
    print m


