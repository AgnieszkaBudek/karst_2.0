#include "network.h"

void::Network::read_setup_file(ifstream& fp_setup){

	while(getline(fp_setup)){


	}



    // parse the setup file
    bool stopper = true;
    while(stopper)
    {
        if(myGetLine(fp_setup, newLine, 2048, &constants, false, &repeatText, readRepeat, repeating, repeatPointer) != EOF)
        {
            int numOfTags = firstWords(newLine, tag, tag2);

            // timestep
            if(sscanf(newLine, "timestep %lf", &dein1) == 1) {
                if(dein1 > 0) {
                    timeStep = dein1;
                    //tutaj beszczelnie wpisze moj krok czasowy
                    //timeStep = pobierz_delta_();
                    std::cerr << "Global setup: timestep = " << timeStep << " ps" << std::endl;
                    p = new ProteinParser(timeStep);
                } else {
                    std::cerr << "Global setup: timestep " << dein1 << " <= 0" << std::endl;
                    exit(-1);
                }
            }
            // simulation length
            else if(sscanf(newLine, "runsteps %ld", &rein1) == 1) {
                if(rein1 > 0) {
                    dauer = rein1;
                    //tutaj beszczelnie wpisze swoj czas symulacji
                    //dauer = pobierz_mstep_();
                    std::cerr << "Global setup: duration = " << dauer << " timesteps = " << timeStep * dauer * 1e-6 << " micro secs" << std::endl;
                } else {
                    std::cerr << "Global setup: simulation length " << rein1 << " <= 0" << std::endl;
                    exit(-1);
                }
            }
            else if(sscanf(newLine, "number %s", helfer) == 1)
            {
                if(strlen(outputFile) > 0) {	// filename is defined
                    // test whether this protein type was already defined
                    list<string>::iterator it;
                    bool isDefined = false;
                    for(it = definedParticles.begin(); it != definedParticles.end(); ++it) {
                        if((*it) == string(helfer)) {
                            isDefined = true;
                            break;
                        }
                    }
                    if(isDefined) {
                        std::cerr << "Global setup: registered protein type >>" << helfer << "<< for output" << std::endl;
                        numbers.push_back(helfer);
                    } else {
                        std::cerr << "Global setup: error registering protein type >>" << helfer << "<< for output (not defined?)" << std::endl;
                        exit(-1);
                    }
                }
                else {
                    std::cerr << "Global setup: no output file defined, cannot register protein type >>" << helfer << "<<" << std::endl;
                    exit(-1);
                }
            }


            // anything that does not fit any of the above tests must be wrong
            else {
                std::cerr << "Global setup: cannot make sense of setup line >>" << newLine << "<<" << std::endl;
                exit(-1);
            }
        }
        // myGetLine returned EOF => We're done, jipii!??
        else
        {
            // see if we're in an included file
            if(fileStack.size() > 0) {
                std::cerr << "Global setup: done with include file\n" << std::endl;
                fclose(fp_setup);
                fp_setup = fileStack.back();
                fileStack.pop_back();
            }
            else {
                std::cerr << "Global setup: reached EOF" << std::endl;
                stopper = false;
                fclose(fp_setup);
            }
        }
    }

}

	return;
}

// Read next line, chop off leading whitespace and ignore comment lines
int myGetLine(FILE *fp_in, char *newLine, const int maxLen, std::map<std::string, std::string> *constants, bool noLowerCase,
              std::list<std::string> *repeatText, const bool readRepeat, const bool repeating, std::list<std::string>::iterator &repPtr)
{
    int c;
    int count = 0;
    bool leadingWhite = true;
    bool ignoreLine = true;
    char myLine [2049];
    char endRepeat[] = "endrepeat";

    // std::cerr << "myGetLine: start reading" << std::endl;

    if(repeating == true)
    {
        // std::cerr << "myGetLine: repeating = true" << std::endl;
        strcpy(myLine, repPtr->c_str());
        repPtr++;
    }
    else
    {
        bool readSilent = true;

        while(readSilent)
        {
            readSilent = false;		// activate switch
            ignoreLine = true;		// have to reset this one

            while(ignoreLine)
            {
                leadingWhite = true;
                count = 0;

                c = fgetc(fp_in);

                // std::cerr << "myGetLine:    count = " << count << ",  c = " << (char)(c) << std::endl;

                while(c != '\n') {
                    if(c == EOF) {
                        return(c);
                    }

                    if((c != ' ' && c != '\t')  || leadingWhite == false) {
                        leadingWhite = false;
                        myLine[count] = c;
                        count++;
                        if(count >= maxLen) {
                            myLine[count-1] = '\0';
                            // return(0);
                            break;
                        }
                    }
                    c = fgetc(fp_in);
                    // std::cerr << "   myGetLine: count = " << count << ",  c = " << (char)(c) << std::endl;
                }

                myLine[count] = '\0';
                if(count > 0 && myLine[0] != '#') {
                    ignoreLine = false;
                }
            }

            // std::cerr << "myGetLine: read myLine = >>" << myLine << "<< from file" << std::endl;


            // save read line into repeat buffer
            if(readRepeat == true) {
                repeatText->push_back(std::string(myLine));

                // std::cerr << "myGetLine: repeatText-buffer:" << std::endl;
                // for(std::list<std::string>::iterator it = repeatText->begin(); it != repeatText->end(); it++) {
                // 	std::cerr << "   " << (*it) << std::endl;
                // }

                // check for "endrepeat"
                int i;
                for(i=0; i<9; i++) {
                    if(tolower(myLine[i]) != endRepeat[i]) {
                        // std::cerr << "myGetLine: comparing >>" << myLine << " to >>endRepeat<<  i = " << i << "  " << tolower(myLine[i]) << " - " << endRepeat[i] << std::endl;
                        break;
                    }
                }
                // std::cerr << "myGetLine: Comparison done, i = "<< i << std::endl;

                if(i<8) {
                    // std::cerr << "  myGetLine: readSilent = true" << std::endl;
                    readSilent = true;
                }
            }
        }
    }


    // std::cerr << "myGetLine: reading done, start processing: myLine = >>" << myLine << "<<" << std::endl;

    // check for replacements
    int i, k, it, itt;
    bool inConst = false;
    char target[2048];

    k = 0;
    i = 0;
    while(myLine[i] != '\0' || inConst)
    {
        char next = myLine[i];

        if(inConst)
        {
            // check whether a replacement is known
            if(constants->count(std::string(target))) {
                for(itt=0; itt<(*constants)[target].size(); itt++) {
                    newLine[k] = (*constants)[target].c_str()[itt];
                    k++;
                }
                inConst = false;
                i--;				// restart from current char
            }
            else  // keep reading target string
            {
                // std::cerr << "myGetLine: no replacement found for constant >>" << target << "<<" << std::endl;

                if(next != ' ' && next != '\t') {
                    target[it] = next;
                    target[it+1] = '\0';
                    it++;
                }
                else {	// replacement not known - exit!
                    std::cerr << "myGetLine: constant >>" << target << "<< is undefined in line >>" << myLine << "<<" << std::endl;
                    exit(-4);
                }
                // std::cerr << "myGetLine: target = >>" << target << "<<" << std::endl;
            }
        }

        else if(next == '$')  {	// start reading target
            inConst = true;
            it = 0;
            target[it] = '\0';
        }
        else {				// just copy
            newLine[k] = next;
            k++;
        }

        i++;
    }
    newLine[k] = '\0';

    // std::cerr << "myGetLine: nearly finished, newLine = >>" << newLine << "<<" << std::endl;

    // convert 1st word to lower case
    i=0;
    while(newLine[i] != ' ' && newLine[i] != '\t') {
        if(!noLowerCase) {
            newLine[i] = tolower(newLine[i]);
        }
        i++;
    }
    // or convert second word to lower for bond types
    if(noLowerCase) {
        // std::cerr << "myGetLine: convert second word tolower(), i = " << i << "  newline[i] = >>" << newLine[i] <<"<<" << std::endl;

        while(newLine[i] == ' ' || newLine[i] == '\t') {
            i++;
        }
        while(newLine[i] != ' ' && newLine[i] != '\t') {
            newLine[i] = tolower(newLine[i]);
            i++;
        }
    }

    // std::cerr << "myGetLine: finished,       newLine = >>" << newLine << "<<" << std::endl;

    return(0);
}
