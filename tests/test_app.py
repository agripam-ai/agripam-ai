import os
import unittest
from pathlib import Path


class AppTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("AGRIPAM_RUN_APP_TESTS") == "1", "set AGRIPAM_RUN_APP_TESTS=1")
    def test_builtin_genome_workflow(self):
        from streamlit.testing.v1 import AppTest

        app = Path(__file__).resolve().parents[1] / "app" / "streamlit_app.py"
        at = AppTest.from_file(str(app), default_timeout=20).run()
        source = next(widget for widget in at.selectbox if widget.label == "Genome source")
        source.select("Built-in demonstration").run()
        submit = next(widget for widget in at.button if widget.label == "Analyze genome")
        submit.click().run(timeout=20)
        self.assertFalse(at.exception)
        metrics = {metric.label: metric.value for metric in at.metric}
        self.assertEqual(metrics["TTC candidates"], "1")


if __name__ == "__main__":
    unittest.main()
