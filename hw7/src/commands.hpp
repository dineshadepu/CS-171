#ifndef COMMANDS_HPP
#define COMMANDS_HPP

#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>

#include "model.hpp"

using namespace std;

static const unsigned int cmd_name_buffer_len = 32;
struct CommandName {
    char name[cmd_name_buffer_len];

    CommandName();  // disabled in definition DO NOT USE
    CommandName(const char* name);

    const bool operator==(const CommandName& rhs) const;
    const bool operator!=(const CommandName& rhs) const;
    const bool operator==(const char* rhs) const;
    const bool operator!=(const char* rhs) const;
};

struct CommandNameHasher {
    unsigned int operator()(const CommandName& cmd) const {
        unsigned int hash = 5381;
        for (unsigned int i = 0; i < cmd_name_buffer_len; i++) {
            hash = ((hash << 5) + hash) + (unsigned int) cmd.name[i];
        }

        return hash;
    }
};

enum CommandType {DEFAULT, STATE};

struct Command {
    CommandType type;

    unsigned int min_expected_argc;
    unsigned int max_expected_argc;
    bool (*action)(int expected_argc, char** argv);

    Command();  // disabled in definition DO NOT USE

    Command(
        CommandType type,
        unsigned int min_expected_argc,
        unsigned int max_expected_argc,
        bool (*action)(int argc, char** argv));
};

struct HelpInfo {
private:
    static const int help_text_buffer_len = 512;

public:
    CommandName cmd_name;
    char help_text[help_text_buffer_len];

    HelpInfo(); // disabled in definition DO NOT USE

    HelpInfo(CommandName cmd_name, const char* help_text);
};

namespace Commands {
    // can't be rigorously checked so be VERY careful not to go over this limit
    // when adding in new commands
    static const int max_cmd_id = 1000;

    // IDs
    // general
    static const int invalid_cmd_id                 = -1;
    static const int help_cmd_id                    = 0;
    static const int stop_cmd_id                    = 1;
    static const int info_cmd_id                    = 2;
    static const int interact_cmd_id                = 3;
    static const int source_cmd_id                  = 4;
    static const int save_cmd_id                    = 5;
    // selection and modifications of Renderables
    static const int deselect_cmd_id                = 99;
    // primitive
    static const int primitive_get_cmd_id           = 100;
    static const int primitive_set_coeff_cmd_id     = 101;
    static const int primitive_set_exponent_cmd_id  = 102;
    static const int primitive_set_patch_cmd_id     = 103;
    static const int primitive_set_color_cmd_id     = 104;
    static const int primitive_set_ambient_cmd_id   = 105;
    static const int primitive_set_reflected_cmd_id = 106;
    static const int primitive_set_refracted_cmd_id = 107;
    static const int primitive_set_gloss_cmd_id     = 108;
    static const int primitive_set_diffuse_cmd_id   = 109;
    static const int primitive_set_specular_cmd_id  = 110;
    // object
    static const int object_get_cmd_id              = 200;
    static const int object_all_translate_cmd_id    = 201;
    static const int object_all_rotate_cmd_id       = 202;
    static const int object_all_scale_cmd_id        = 203;
    static const int object_add_object_cmd_id       = 204;
    static const int object_add_primitive_cmd_id    = 205;
    static const int object_set_cursor_cmd_id       = 206;
    static const int object_cursor_translate_cmd_id = 207;
    static const int object_cursor_rotate_cmd_id    = 208;
    static const int object_cursor_scale_cmd_id     = 209;

    // name to ID definitions
    static const unordered_map<CommandName, int, CommandNameHasher>
        name_to_id({
            // general
            {CommandName("help"),           help_cmd_id},
            {CommandName("quit"),           stop_cmd_id},
            {CommandName("exit"),           stop_cmd_id},
            {CommandName("stop"),           stop_cmd_id},
            {CommandName("info"),           info_cmd_id},
            {CommandName("ls"),             info_cmd_id},
            {CommandName("interact"),       interact_cmd_id},
            {CommandName("source"),         source_cmd_id},
            {CommandName("load"),           source_cmd_id},
            {CommandName("save"),           save_cmd_id},
            // selection and modifications of Renderables
            {CommandName("deselect"),       deselect_cmd_id},
            // primitive
            {CommandName("Primitive"),      primitive_get_cmd_id},
            {CommandName("primitive"),      primitive_get_cmd_id},
            {CommandName("prm"),            primitive_get_cmd_id},
            {CommandName("setCoeff"),       primitive_set_coeff_cmd_id},
            {CommandName("Coeff"),          primitive_set_coeff_cmd_id},
            {CommandName("coeff"),          primitive_set_coeff_cmd_id},
            {CommandName("setExponent"),    primitive_set_exponent_cmd_id},
            {CommandName("exponent"),       primitive_set_exponent_cmd_id},
            {CommandName("Exponent"),       primitive_set_exponent_cmd_id},
            {CommandName("exp"),            primitive_set_exponent_cmd_id},
            {CommandName("setPatch"),       primitive_set_patch_cmd_id},
            {CommandName("Patch"),          primitive_set_patch_cmd_id},
            {CommandName("patch"),          primitive_set_patch_cmd_id},
            {CommandName("setColor"),       primitive_set_color_cmd_id},
            {CommandName("color"),          primitive_set_color_cmd_id},
            {CommandName("surface"),        primitive_set_color_cmd_id},
            {CommandName("setAmbient"),     primitive_set_ambient_cmd_id},
            {CommandName("ambient"),        primitive_set_ambient_cmd_id},
            {CommandName("Ambient"),        primitive_set_ambient_cmd_id},
            {CommandName("setReflected"),   primitive_set_reflected_cmd_id},
            {CommandName("Reflected"),      primitive_set_reflected_cmd_id},
            {CommandName("reflected"),      primitive_set_reflected_cmd_id},
            {CommandName("setRefracted"),   primitive_set_refracted_cmd_id},
            {CommandName("Refracted"),      primitive_set_refracted_cmd_id},
            {CommandName("refracted"),      primitive_set_refracted_cmd_id},
            {CommandName("setGloss"),       primitive_set_gloss_cmd_id},
            {CommandName("gloss"),          primitive_set_gloss_cmd_id},
            {CommandName("setDiffuse"),     primitive_set_diffuse_cmd_id},
            {CommandName("diffuse"),        primitive_set_diffuse_cmd_id},
            {CommandName("setSpecular"),    primitive_set_specular_cmd_id},
            {CommandName("specular"),       primitive_set_specular_cmd_id},
            // object
            {CommandName("Object"),         object_get_cmd_id},
            {CommandName("object"),         object_get_cmd_id},
            {CommandName("obj"),            object_get_cmd_id},
            {CommandName("translateAll"),   object_all_translate_cmd_id},
            {CommandName("xall"),           object_all_translate_cmd_id},
            {CommandName("rotateAll"),      object_all_rotate_cmd_id},
            {CommandName("rall"),           object_all_rotate_cmd_id},
            {CommandName("scaleAll"),       object_all_scale_cmd_id},
            {CommandName("sall"),           object_all_scale_cmd_id},
            {CommandName("addObject"),      object_add_object_cmd_id},
            {CommandName("ao"),             object_add_object_cmd_id},
            {CommandName("addPrimitive"),   object_add_primitive_cmd_id},
            {CommandName("ap"),             object_add_primitive_cmd_id},
            {CommandName("setCursor"),      object_set_cursor_cmd_id},
            {CommandName("cursor"),         object_set_cursor_cmd_id},
            {CommandName("cur"),            object_set_cursor_cmd_id},
            {CommandName("translate"),      object_cursor_translate_cmd_id},
            {CommandName("xs"),             object_cursor_translate_cmd_id},
            {CommandName("rotate"),         object_cursor_rotate_cmd_id},
            {CommandName("rs"),             object_cursor_rotate_cmd_id},
            {CommandName("scale"),          object_cursor_scale_cmd_id},
            {CommandName("ss"),             object_cursor_scale_cmd_id}
        });

    // command callback functions
    // general
    bool help(int argc, char** argv);
    bool stop(int argc, char** argv);
    bool info(int argc, char** argv);
    bool interact(int argc, char** argv);
    bool source(int argc, char** argv);
    bool save(int argc, char** argv);
    // selection and modifications of Renderables
    bool deselect(int argc, char** argv);
    // primitive
    bool getPrimitive(int argc, char** argv);
    bool prmSetCoeff(int argc, char** argv);
    bool prmSetExponent(int argc, char** argv);
    bool prmSetPatch(int argc, char** argv);
    bool prmSetColor(int argc, char** argv);
    bool prmSetAmbient(int argc, char** argv);
    bool prmSetReflected(int argc, char** argv);
    bool prmSetRefracted(int argc, char** argv);
    bool prmSetGloss(int argc, char** argv);
    bool prmSetDiffuse(int argc, char** argv);
    bool prmSetSpecular(int argc, char** argv);
    // object
    bool getObject(int argc, char** argv);
    bool objAllTranslate(int argc, char** argv);
    bool objAllRotate(int argc, char** argv);
    bool objAllScale(int argc, char** argv);
    bool objAddObject(int argc, char** argv);
    bool objAddPrimitive(int argc, char** argv);
    bool objSetCursor(int argc, char** argv);
    bool objCursorTranslate(int argc, char** argv);
    bool objCursorRotate(int argc, char** argv);
    bool objCursorScale(int argc, char** argv);

    // ID to command definitions
    static unordered_map<int, Command> id_to_cmd({
        // general
        {invalid_cmd_id,                 Command(DEFAULT, 0, 0, NULL)},
        {help_cmd_id,                    Command(DEFAULT, 1, (unsigned int) -1, &help)},
        {stop_cmd_id,                    Command(DEFAULT, 1, 1, &stop)},
        {info_cmd_id,                    Command(DEFAULT, 1, (unsigned int) -1, &info)},
        {interact_cmd_id,                Command(DEFAULT, 1, 1, &interact)},
        {source_cmd_id,                  Command(DEFAULT, 1, 2, &source)},
        {save_cmd_id,                    Command(DEFAULT, 1, 2, &save)},
        // selection and modifications of Renderables
        {deselect_cmd_id,                Command(DEFAULT, 1, 1, &deselect)},
        // primitive
        {primitive_get_cmd_id,           Command(STATE, 2, 2, &getPrimitive)},
        {primitive_set_coeff_cmd_id,     Command(DEFAULT, 4, 4, &prmSetCoeff)},
        {primitive_set_exponent_cmd_id,  Command(DEFAULT, 3, 3, &prmSetExponent)},
        {primitive_set_patch_cmd_id,     Command(DEFAULT, 3, 3, &prmSetPatch)},
        {primitive_set_color_cmd_id,     Command(DEFAULT, 4, 4, &prmSetColor)},
        {primitive_set_ambient_cmd_id,   Command(DEFAULT, 2, 2, &prmSetAmbient)},
        {primitive_set_reflected_cmd_id, Command(DEFAULT, 2, 2, &prmSetReflected)},
        {primitive_set_refracted_cmd_id, Command(DEFAULT, 2, 2, &prmSetRefracted)},
        {primitive_set_gloss_cmd_id,     Command(DEFAULT, 2, 2, &prmSetGloss)},
        {primitive_set_diffuse_cmd_id,   Command(DEFAULT, 2, 2, &prmSetDiffuse)},
        {primitive_set_specular_cmd_id,  Command(DEFAULT, 2, 2, &prmSetSpecular)},
        // object
        {object_get_cmd_id,              Command(STATE, 2, 2, &getObject)},
        {object_all_translate_cmd_id,    Command(DEFAULT, 4, 4, &objAllTranslate)},
        {object_all_rotate_cmd_id,       Command(DEFAULT, 5, 5, &objAllRotate)},
        {object_all_scale_cmd_id,        Command(DEFAULT, 4, 4, &objAllScale)},
        {object_add_object_cmd_id,       Command(DEFAULT, 2, 3, &objAddObject)},
        {object_add_primitive_cmd_id,    Command(DEFAULT, 2, 3, &objAddPrimitive)},
        {object_set_cursor_cmd_id,       Command(DEFAULT, 2, 2, &objSetCursor)},
        {object_cursor_translate_cmd_id, Command(DEFAULT, 4, 4, &objCursorTranslate)},
        {object_cursor_rotate_cmd_id,    Command(DEFAULT, 5, 5, &objCursorRotate)},
        {object_cursor_scale_cmd_id,     Command(DEFAULT, 4, 4, &objCursorScale)}
    });

    static unordered_map<int, HelpInfo*> help_infos({
        // general
        {stop_cmd_id, new HelpInfo(CommandName("stop"), "\
behavior:           exits the program. does not save progress. same behavior as\n\
                    ctrl+c\n\
expected arguments: [NONE]\n")},

        {info_cmd_id, new HelpInfo(CommandName("info"), "\
behavior:           displays information on active Renderables. if no Renderable\n\
                    is specified scene information is displayed instead.\n\
expected arguments: [Renderable Name(s)]\n\n\
NOTE: passing 'all' as an argument will display scene information.\n\
passing 'selected' as an argument will display information on currently selected\n\
Renderable\n")},

        {interact_cmd_id, new HelpInfo(CommandName("interact"), "\
behavior:           temporarily disables commandline and enables mouse and\n\
                    interaction with the display. press q or esc at the display\n\
                    to return to commandline.\n\
expected arguments: [NONE]\n")},

        {source_cmd_id, new HelpInfo(CommandName("source"), "\
behavior:           loads a scene stored in a .scn file. if no file is\n\
                    specified, searches for data/quicksave.scn instead. if no\n\
                    file can be found and loaded, no action is taken.\n\
expected arguments: [filename]\n\n\
NOTE: loading a file will clear all current scene information. make sure to save\n\
any progress before loading a scene\n")},

        {save_cmd_id, new HelpInfo(CommandName("save"), "\
behavior:           saves the current scene into a .scn file in the data\n\
                    directory. if no filename is specified, the current scene is\n\
                    saved into data/quicksave.scn\n\
expected arguments: [filename]\n\n\
NOTE: saving to a pre-existing file will overwrite the data\n")},

        // selection and modifications of Renderables
        {deselect_cmd_id, new HelpInfo(CommandName("deselect"), "\
behavior:           removes any Renderable selection. will not affect the scene\n\
expected arguments: [NONE]\n")},

        // primitive
        {primitive_get_cmd_id, new HelpInfo(CommandName("Primitive"), "\
behavior:           selects a Primitive generating one if one with the given\n\
                    name does not exist yet\n\
expected arguments: [name]\n\n\
NOTE: all Renderables share the same name pool which means that primitives and\n\
objects cannot have the same names\n")},

        {primitive_set_coeff_cmd_id, new HelpInfo(CommandName("setCoeff"), "\
behavior:           sets the scaling coefficients of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [x] [y] [z] (all args must be non negative numeric values)\n")},

        {primitive_set_exponent_cmd_id, new HelpInfo(CommandName("setExponent"), "\
behavior:           sets the exponents of the currently selected primitive. must\n\
                    have a primitive selected\n\
expected arguments: [e] [n] (all args must be non negative numeric values\n")},

        {primitive_set_patch_cmd_id, new HelpInfo(CommandName("setPatch"), "\
behavior:           sets the patch count of the currently selected primitive.\n\
                    must have a primitive selected\n\
expected arguments: [longitudinal patch count] [latitudinal patch count] (all\n\
                    args must be non negative integral values)\n")},

        {primitive_set_color_cmd_id, new HelpInfo(CommandName("setColor"), "\
behavior:           sets the color of the currently selected primitive. must\n\
                    a primitive selected\n\
expected arguments: [r] [g] [b] (all args must be numeric values between 0.0 and\n\
                    1.0 (inclusive))\n")},

        {primitive_set_ambient_cmd_id, new HelpInfo(CommandName("setAmbient"), "\
behavior:           sets the ambient color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [a] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        {primitive_set_reflected_cmd_id, new HelpInfo(CommandName("setReflected"), "\
behavior:           sets the reflected color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [r] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        {primitive_set_refracted_cmd_id, new HelpInfo(CommandName("setRefracted"), "\
behavior:           sets the refracted color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [r] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        {primitive_set_gloss_cmd_id, new HelpInfo(CommandName("setGloss"), "\
behavior:           sets the gloss color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [g] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        {primitive_set_diffuse_cmd_id, new HelpInfo(CommandName("setDiffuse"), "\
behavior:           sets the diffuse color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [d] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        {primitive_set_specular_cmd_id, new HelpInfo(CommandName("setSpecular"), "\
behavior:           sets the specular color scale of the currently selected\n\
                    primitive. must have a primitive selected\n\
expected arguments: [s] (a single numeric value between 0.0 and 1.0 (inclusive)\n")},

        // object
        {object_get_cmd_id, new HelpInfo(CommandName("Object"), "\
behavior:           selects an Object generating one if one with the given\n\
                    name does not exist yet\n\
expected arguments: [name]\n\n\
NOTE: all Renderables share the same name pool which means that primitives and\n\
objects cannot have the same names\n")},

        {object_all_translate_cmd_id, new HelpInfo(CommandName("translateAll"), "\
behavior:           applies a translation to the currently selected object. must\n\
                    have an Object selected\n\
expected arguments: [x] [y] [z] (all args must be numeric values)\n")},

        {object_all_rotate_cmd_id, new HelpInfo(CommandName("rotateAll"), "\
behavior:           applies a rotation to the currently selected object. must\n\
                    have an Object selected\n\
expected arguments: [x] [y] [z] [theta] (all args must be numeric values and\n\
                    theta is in degrees)\n\n\
NOTE: the axis given by [x] [y] [z] does not have to be normalized\n")},

        {object_all_scale_cmd_id, new HelpInfo(CommandName("scaleAll"), "\
behavior:           applies a scaling transformation to the currently selected\n\
                    object. must have an Object selected\n\
expected arguments: [x] [y] [z] (all args must be numeric values)\n")},

        {object_add_object_cmd_id, new HelpInfo(CommandName("addObject"), "\
behavior:           adds a child Object to the currently selected Object. must\n\
                    have an Object selected\n\
expected arguments: [name]\n")},

        {object_add_primitive_cmd_id, new HelpInfo(CommandName("addPrimitive"), "\
behavior:           adds a child Primitive to the currently selected Object.\n\
                    must have an Object selected\n\
expected arguments: [name]\n")},

        {object_set_cursor_cmd_id, new HelpInfo(CommandName("setCursor"), "\
behavior:           sets child cursor of the currently selected Object. must\n\
                    have an Object selected\n\
expected arguments: [name]\n")},

        {object_cursor_translate_cmd_id, new HelpInfo(CommandName("translate"), "\
behavior:           applies a translation to the child pointed to by the cursor\n\
                    of the currently selected Object. must have an Object with\n\
                    at least one child selected\n\
expected arguments: [x] [y] [z]\n")},

        {object_cursor_rotate_cmd_id, new HelpInfo(CommandName("rotate"), "\
behavior:           applies a rotation to the child pointed to by the cursor of\n\
                    the currently selected Object. must have an Object with at\n\
                    least one child selected\n\
expected arguments: [x] [y] [z] [theta] (all args must be numeric values and\n\
                    theta is in degrees)\n\n\
NOTE: the axis given by [x] [y] [z] does not have to be normalized\n")},

        {object_cursor_scale_cmd_id, new HelpInfo(CommandName("scale"), "\
behavior:           applies a scaling transformation to the child pointed to by\n\
                    the cursor of the currently selected Object. must have an\n\
                    Object with at least one child selected\n\
expected arguments: [x] [y] [z]\n")},
    });

    void init();
    void help_init();
};

#endif