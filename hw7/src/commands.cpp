#include "commands.hpp"
#include "command_line.hpp"

/******************************** CommandName *********************************/

CommandName::CommandName() {
    // disable
    fprintf(stderr, "ERROR CommandName is not allowed a default construction\n");
    exit(1);
}
CommandName::CommandName(const char* name) : name() {
    strcpy(this->name, name);
}

const bool CommandName::operator==(const CommandName& rhs) const {
    return strcmp(this->name, rhs.name) == 0;
}
const bool CommandName::operator!=(const CommandName& rhs) const {
    return strcmp(this->name, rhs.name) != 0;
}
const bool CommandName::operator==(const char* rhs) const {
    return strcmp(this->name, rhs) == 0;
}
const bool CommandName::operator!=(const char* rhs) const {
    return strcmp(this->name, rhs) != 0;
}

/********************************** Command ***********************************/

Command::Command() {
    // disable
    fprintf(stderr, "ERROR Command is not allowed a default construction\n");
    exit(1);
}

Command::Command(
    CommandType type,
    unsigned int min_expected_argc,
    unsigned int max_expected_argc,
    bool (*action)(int argc, char** argv)) :
    type(type),
    min_expected_argc(min_expected_argc),
    max_expected_argc(max_expected_argc),
    action(action)
{

}

/********************************* HelpInfo ***********************************/

HelpInfo::HelpInfo() {
    // disable
    fprintf(stderr, "ERROR HelpInfo is not allowed a default construction\n");
    exit(1);
}

HelpInfo::HelpInfo(CommandName cmd_name, const char* help_text) :
    cmd_name(cmd_name),
    help_text()
{
    strcpy(this->help_text, help_text);
}













/***************************** Command Callbacks ******************************/

// general
bool Commands::help(int argc, char** argv) {
    assert(argc > 0);
    
    if (argc == 1) {
        // print all possible commands
        printf("all possible commands:\n");
        for (int i = 0; i <= Commands::max_cmd_id; i++) {
            if (Commands::help_infos.find(i) != Commands::help_infos.end()) {
                const HelpInfo* help_info = Commands::help_infos[i];
                printf("\t%s", help_info->cmd_name.name);
                int cmd_len = strlen(help_info->cmd_name.name);
                for (int i = 0; i < 15 - cmd_len; i++) {
                    printf(" ");
                }
                vector<CommandName> aliases;
                for (const auto& cmd : Commands::name_to_id) {
                    if (cmd.second == i && cmd.first != help_info->cmd_name) {
                        aliases.push_back(cmd.first);
                    }
                }
                if (aliases.size() != 0) {
                    printf("(aliases:");
                    for (auto& alias : aliases) {
                        printf(" %s", alias.name);
                    }
                    printf(")");
                }
                printf("\n");
            }
        }
        printf("for detail(s) enter 'help [command name(s) or alias(es)]'\n");
    } else {
        for (int i = 1; i < argc; i++) {
            const CommandName& cmd_name(argv[i]);
            if (Commands::name_to_id.find(cmd_name) ==
                Commands::name_to_id.end())
            {
                printf("%s is not a valid command or alias\n", argv[i]);
            } else {
                int cmd_id = Commands::name_to_id.at(cmd_name);
                printf("%s", Commands::help_infos[cmd_id]->help_text);
            }
        }
    }
    return false;
}
bool Commands::stop(int argc, char** argv) {
    exit(0);
}
bool Commands::info(int argc, char** argv) {
    assert(argc > 0);
    const Line* cur_state = CommandLine::getState();

    if (argc == 1) {
        // print scene information
        printSceneInfo(0);

        // print current selection
        if (cur_state) {
            int cmd_id = cur_state->toCommandID();
            switch (cmd_id) {
                case Commands::primitive_get_cmd_id:
                    printf("currently selected: PRM %s\n",
                        cur_state->tokens[1]);
                    break;
                case Commands::object_get_cmd_id:
                    printf("currently selected: OBJ %s\n",
                        cur_state->tokens[1]);
                    break;
                default:
                    fprintf(stderr, "ERROR Commands:info invalid state CommandID %d from current state\n",
                        cmd_id);
                    exit(1);
            }
        } else {
            printf("no Renderable currently selected\n");
        }
    } else {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "all") == 0) {
                printSceneInfo(0);
            } else if (strcmp(argv[i], "selected") == 0) {
                // print current selection
                if (cur_state) {
                    int cmd_id = cur_state->toCommandID();
                    switch (cmd_id) {
                        case Commands::primitive_get_cmd_id:
                            printf("currently selected: PRM %s\n",
                                cur_state->tokens[1]);
                            break;
                        case Commands::object_get_cmd_id:
                            printf("currently selected: OBJ %s\n",
                                cur_state->tokens[1]);
                            break;
                        default:
                            fprintf(stderr, "ERROR Commands:info invalid state CommandID %d from current state\n",
                                cmd_id);
                            exit(1);
                    }
                    printInfo(Renderable::get(cur_state->tokens[1]), 1);
                } else {
                    printf("no Renderable currently selected\n");
                }
            } else {
                if (Renderable::exists(argv[i])) {
                    printf("%s:\n", argv[i]);
                    printInfo(Renderable::get(argv[i]), 1);
                } else {
                    fprintf(stderr, "printInfo ERROR Renderable with name %s does not exist",
                        argv[i]);
                }
            }
            printf("\n");
        }
    }
    return false;
}
bool Commands::interact(int argc, char** argv) {
    assert(argc == 1);
    CommandLine::pause();
    return false;
}
bool Commands::source(int argc, char** argv) {
    assert(argc == 1 || argc == 2);

    bool good_response = false;
    bool approved = false;
    printf("current scene and history will be cleared if you load a file. continue (yes/no)? ");
    char input_buffer[64];
    while (!good_response) {
        cin.getline(input_buffer, 64);
        if (strcmp(input_buffer, "yes") == 0) {
            approved = true;
            good_response = true;
        } else if (strcmp(input_buffer, "no") == 0) {
            good_response = true;
        } else {
            printf("invalid response. please enter 'yes' or 'no' ");
        }
    }

    if (approved) {
        Renderable::clear();
        CommandLine::clearState();
        CommandLine::clearHistory();

        ifstream source_file;

        if (argc == 1) {
            source_file.open(quicksave_dir, ifstream::in);
        } else {
            source_file.open(argv[1], ifstream::in);
        }
        if (source_file.is_open()) {
            while (!source_file.eof()) {
                CommandLine::readLine(source_file);
            }

            source_file.close();
        } else {
            fprintf(stderr, "ERROR couldn't open file %s\n", argv[1]);
        }
    } else {
        printf("source command aborted.\n");
    }
    return false;
}
bool Commands::save(int argc, char** argv) {
    assert(argc == 1 || argc == 2);

    ofstream savefile;

    if (argc == 1) {
        savefile.open(quicksave_dir);
    } else {
        assert(strlen(argv[1]) < 256);
        char filename_buffer[256];
        strcpy(filename_buffer, argv[1]);
        savefile.open(filename_buffer);
    }

    for (const auto& line : CommandLine::getHistory()) {
        savefile << line->untokenized_argv << endl;
    }
    savefile.close();
    return false;
}
// selection and modifications of Renderables
bool Commands::deselect(int argc, char** argv) {
    assert(argc == 1);
    CommandLine::clearState();
    return true;
}
// primitive
bool Commands::getPrimitive(int argc, char** argv) {
    Renderable* ren = Renderable::get(argv[1]);
    if (ren) {
        if (ren->getType() != PRM) {
            fprintf(stderr, "Commands::getPrimitive ERROR Renderable with name %s already exists and has type %s\n",
                argv[1], toCstr(ren->getType()));
            CommandLine::clearState();
            return false;
        }
    } else {
        Renderable::create(PRM, argv[1]);
    }
    return true;
}
bool Commands::prmSetCoeff(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            if (x > 0.0 && y > 0.0 && z > 0.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setCoeff(x, y, z);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 3 non negative numeric values\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetExponent(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float Exp0 = atof(argv[1]);
            float Exp1 = atof(argv[2]);
            if (Exp0 > 0.0 && Exp1 > 0.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setExponents(Exp0, Exp1);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 2 non negative numeric values\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetPatch(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            int patch_x = atoi(argv[1]);
            int patch_y = atoi(argv[2]);
            float patch_x_f = atof(argv[1]);
            float patch_y_f = atof(argv[2]);
            if (patch_x_f > 0.0 && patch_y_f > 0.0 &&
                patch_x_f - patch_x == 0 && patch_y_f - patch_y == 0)
            {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setPatch(patch_x, patch_y);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 2 non negative integers\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetColor(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float r = atof(argv[1]);
            float g = atof(argv[2]);
            float b = atof(argv[3]);
            if (r >= 0.0 && g >= 0.0 && b >= 0.0 &&
                r <= 1.0 && g <= 1.0 && b <= 1.0)
            {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setColor(r, g, b);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 3 numerical values between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetAmbient(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float ambient = atof(argv[1]);
            if (ambient >= 0.0 && ambient <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setAmbient(ambient);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetReflected(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float reflected = atof(argv[1]);
            if (reflected >= 0.0 && reflected <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setReflected(reflected);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetRefracted(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float refracted = atof(argv[1]);
            if (refracted >= 0.0 && refracted <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setRefracted(refracted);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetGloss(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float gloss = atof(argv[1]);
            if (gloss >= 0.0 && gloss <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setGloss(gloss);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetDiffuse(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float diffuse = atof(argv[1]);
            if (diffuse >= 0.0 && diffuse <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setDiffuse(diffuse);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::prmSetSpecular(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::primitive_get_cmd_id) {
            float specular = atof(argv[1]);
            if (specular >= 0.0 && specular <= 1.0) {
                Primitive* prm = dynamic_cast<Primitive*>(
                    Renderable::get(cur_state->tokens[1]));
                prm->setSpecular(specular);
                return true;
            } else {
                fprintf(stderr, "ERROR input values for %s must be 1 numerical value between 0.0 and 1.0 (inclusive)\n",
                    argv[0]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Primitive selected\n",
            argv[0]);
    }
    return false;
}
// object
bool Commands::getObject(int argc, char** argv) {
    Renderable* ren = Renderable::get(argv[1]);
    if (ren) {
        if (ren->getType() != OBJ) {
            fprintf(stderr, "Commands::getObject ERROR Renderable with name %s already exists and has type %s\n",
                argv[1], toCstr(ren->getType()));
            CommandLine::clearState();
            return false;
        }
    } else {
        Renderable::create(OBJ, argv[1]);
    }
    return true;
}
bool Commands::objAllTranslate(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->overallTranslate(x, y, z);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objAllRotate(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            float theta = atof(argv[4]) * M_PI / 180.0;
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->overallRotate(x, y, z, theta);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objAllScale(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->overallScale(x, y, z);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objAddObject(int argc, char** argv) {
    assert(argc == 2 || argc == 3);
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            Renderable* new_child = Renderable::get(argv[1]);
            if (new_child) {
                if (new_child->getType() == OBJ) {
                    Object* obj = dynamic_cast<Object*>(
                        Renderable::get(cur_state->tokens[1]));
                    if (obj->aliasExists(argv[argc - 1])) {
                        fprintf(stderr, "objAddObject ERROR child with alias %s already exists\n",
                            argv[argc - 1]);
                    } else {
                        obj->addChild(argv[1], argv[argc - 1]);
                        return true;
                    }
                } else {
                    fprintf(stderr, "objAddObject ERROR Renderable with name %s has type %s\n",
                        argv[1], toCstr(new_child->getType()));
                }
            } else {
                fprintf(stderr, "objAddObject ERROR Renderable with name %s does not exist\n",
                    argv[1]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objAddPrimitive(int argc, char** argv) {
    assert(argc == 2 || argc == 3);
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            Renderable* new_child = Renderable::get(argv[1]);
            if (new_child) {
                if (new_child->getType() == PRM) {
                    Object* obj = dynamic_cast<Object*>(
                        Renderable::get(cur_state->tokens[1]));
                    if (obj->aliasExists(argv[argc - 1])) {
                        fprintf(stderr, "objAddPrimitive ERROR child with alias %s already exists\n",
                            argv[argc - 1]);
                    } else {
                        obj->addChild(argv[1], argv[argc - 1]);
                        return true;
                    }
                } else {
                    fprintf(stderr, "objAddPrimitive ERROR Renderable with name %s has type %s\n",
                        argv[1], toCstr(new_child->getType()));
                }
            } else {
                fprintf(stderr, "objAddPrimitive ERROR Renderable with name %s does not exist\n",
                    argv[1]);
            }
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objSetCursor(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->setCursor(argv[1]);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objCursorTranslate(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->cursorTranslate(x, y, z);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objCursorRotate(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            float theta = atof(argv[4]) * M_PI / 180.0;
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->cursorRotate(x, y, z, theta);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}
bool Commands::objCursorScale(int argc, char** argv) {
    const Line* cur_state = CommandLine::getState();

    // if state is blank we have nothing selected
    // print error message and skip
    if (cur_state) {
        if (cur_state->toCommandID() == Commands::object_get_cmd_id) {
            float x = atof(argv[1]);
            float y = atof(argv[2]);
            float z = atof(argv[3]);
            // WARNING figure out a way to check that input are numeric values
            Object* obj = dynamic_cast<Object*>(
                Renderable::get(cur_state->tokens[1]));
            obj->cursorScale(x, y, z);
            return true;
        } else {
            fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
                argv[0]);
        }
    } else {
        fprintf(stderr, "ERROR %s requires that you have a Object selected\n",
            argv[0]);
    }
    return false;
}

/************************************** init **********************************/

void Commands::init() {
    help_init();
}

void Commands::help_init() {
}