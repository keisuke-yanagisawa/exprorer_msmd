import tempfile
import jinja2
import os

from .variable import Path


def render_str(template_name: str,
               **kwargs) -> str:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__) + "/.template"))
    template = env.get_template(template_name)
    return template.render(kwargs)


def render_file(template_name: str,
                **kwargs) -> Path:
    out = Path(tempfile.mkstemp(suffix=".in")[1])
    with open(out.path, "w") as fout:
        fout.write(render_str(template_name, **kwargs))
    return out
